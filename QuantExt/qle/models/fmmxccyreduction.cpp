/*
 Copyright (C) 2026 Ed Herranz
 All rights reserved.

 This file is part of ORE, a free-software/open-source library
 for transparent pricing and risk analysis - http://opensourcerisk.org

 ORE is free software: you can redistribute it and/or modify it
 under the terms of the Modified BSD License.  You should have received a
 copy of the license along with this program.
 The license is also available online at <http://opensourcerisk.org>

 This program is distributed on the basis that it will form a useful
 contribution to risk analytics and model standardisation, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the license for more details.
*/

#include <qle/cashflows/overnightindexedcoupon.hpp>
#include <qle/instruments/rebatedexercise.hpp>
#include <qle/models/fmmxccyreduction.hpp>
#include <qle/pricingengines/fmmlsmmultilegoptionengine.hpp>

#include <ql/cashflows/coupon.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/cashflows/overnightindexedcoupon.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/settings.hpp>

#include <cmath>
#include <map>
#include <sstream>

namespace QuantExt {

std::set<Date> fmmXccyContractualDates(const MultiLegOption::arguments& args, const Date& today) {
    std::set<Date> dates;
    for (const auto& leg : args.legs) {
        for (const auto& cf : leg) {
            if (cf->date() <= today)
                continue;
            dates.insert(cf->date());
            if (auto cpn = QuantLib::ext::dynamic_pointer_cast<Coupon>(cf)) {
                if (cpn->accrualStartDate() > today)
                    dates.insert(cpn->accrualStartDate());
                dates.insert(cpn->accrualEndDate());
            }
        }
    }
    if (args.exercise) {
        for (const Date& d : args.exercise->dates())
            if (d > today)
                dates.insert(d);
        for (const Date& d : args.settlementDates)
            if (d > today)
                dates.insert(d);
    }
    return dates;
}

namespace {

// an A floating coupon after telescoping: its accrual start (the +N flow), its payment (the -N
// flow and the spread), all per Hagan eq. 2.2 / 1.4
struct TelescopedCoupon {
    Size startIdx = 0, payIdx = 0;
    Date start, pay;
    Real nominal = 0.0, spread = 0.0, accrual = 0.0, sign = 0.0;
};

struct BFixedCoupon {
    Size startIdx = 0, payIdx = 0;
    Date pay;
    Real nominal = 0.0, accrual = 0.0, signedAmount = 0.0, sign = 0.0;
};

template <class OnCoupon>
void checkOvernightCoupon(const OnCoupon& on, const Date& today, const Handle<YieldTermStructure>& discountA,
                          const std::string& what) {
    QL_REQUIRE(on.accrualStartDate() >= today,
               "fmmMapXccyReduced: " << what << ": coupon accruing at the valuation date is not supported");
    QL_REQUIRE(std::fabs(on.gearing() - 1.0) < 1e-12, "fmmMapXccyReduced: " << what << ": gearing must be 1");
    // the telescoping identity holds when the index projects on the discount curve of its currency
    auto ibor = QuantLib::ext::dynamic_pointer_cast<IborIndex>(on.index());
    QL_REQUIRE(ibor, "fmmMapXccyReduced: " << what << ": floating index is not an ibor / overnight index");
    const auto fwd = ibor->forwardingTermStructure();
    QL_REQUIRE(!fwd.empty(), "fmmMapXccyReduced: " << what << ": floating index without a forwarding curve");
    const Real pFwd = fwd->discount(on.accrualStartDate()) / fwd->discount(on.accrualEndDate());
    const Real pDisc = discountA->discount(on.accrualStartDate()) / discountA->discount(on.accrualEndDate());
    QL_REQUIRE(std::fabs(pFwd / pDisc - 1.0) < 1e-8,
               "fmmMapXccyReduced: " << what << ": the floating index's forwarding curve differs from the discount "
                                        "curve of its currency (period " << on.accrualStartDate() << " - "
                                     << on.accrualEndDate() << ": " << pFwd << " vs " << pDisc
                                     << "); the telescoping identity (Hagan eq. 2.2) needs the discount curve's "
                                        "own overnight index (a deterministic basis must be folded into the spread)");
}

} // namespace

FmmCallableInstrument fmmMapXccyReduced(const MultiLegOption::arguments& args, const FmmGrid& grid,
                                        const Natural toleranceDays, const Date& today,
                                        const FmmXccyReductionInputs& inputs, FmmXccyReductionDiagnostics& diag,
                                        std::vector<Date>& usedSettle) {
    QL_REQUIRE(args.legs.size() == args.payer.size() && args.legs.size() == args.currency.size(),
               "fmmMapXccyReduced: inconsistent legs / payer / currency vectors");
    QL_REQUIRE(args.settlementType == Settlement::Physical,
               "fmmMapXccyReduced: physical settlement only (cash settlement is not supported)");
    QL_REQUIRE(!args.midCouponExercise, "fmmMapXccyReduced: mid-coupon exercise is not supported");
    QL_REQUIRE(inputs.ccyA != inputs.ccyB, "fmmMapXccyReduced: currencies A and B must differ");
    QL_REQUIRE(!inputs.discountA.empty() && !inputs.discountB.empty() && !inputs.fxSpotBperA.empty(),
               "fmmMapXccyReduced: discount curves and the FX spot are required");
    const Size M = grid.numberOfRates();
    diag = FmmXccyReductionDiagnostics();
    diag.fxSpotBperA = inputs.fxSpotBperA->value();
    QL_REQUIRE(diag.fxSpotBperA > 0.0, "fmmMapXccyReduced: non-positive FX spot");

    // F^{A/B}(0; T) = F^{A/B}(0) Z^A(0; T) / Z^B(0; T) (eq. 1.7b) with the pricing curves
    auto fwdFx = [&](const Date& d) {
        return diag.fxSpotBperA * inputs.discountA->discount(d) / inputs.discountB->discount(d);
    };

    FmmCallableInstrument inst;
    inst.style = FmmCallableInstrument::Style::Enter;
    inst.fixedFlows.assign(M + 1, 0.0);
    inst.floatWeights.assign(M + 1, 0.0);
    FmmLegMapper mapper(grid, toleranceDays, today);

    std::vector<std::pair<Size, Size>> couponSpans; // (accrual start index, pay index) of every future coupon
    std::vector<TelescopedCoupon> aFloats;
    std::vector<BFixedCoupon> bFixed;
    std::vector<Real> aStartFlowsB(M + 1, 0.0); // +N at the accrual start of an A floating coupon, in B
    Real lastAPay = 0.0, lastBPay = 0.0;

    for (Size l = 0; l < args.legs.size(); ++l) {
        const Real sign = args.payer[l] ? -1.0 : 1.0;
        std::ostringstream what;
        what << "leg " << l;
        if (args.currency[l] == inputs.ccyB) {
            mapper.addLeg(args.legs[l], sign, inst, what.str());
            for (const auto& cf : args.legs[l]) {
                if (cf->date() <= today)
                    continue;
                if (auto fx = QuantLib::ext::dynamic_pointer_cast<FixedRateCoupon>(cf)) {
                    BFixedCoupon c;
                    c.payIdx = grid.index(fx->date(), toleranceDays, what.str() + " fixed pay date");
                    c.startIdx = fx->accrualStartDate() > today
                                     ? grid.index(fx->accrualStartDate(), toleranceDays, what.str() + " accrual start")
                                     : 0;
                    c.pay = fx->date();
                    c.nominal = fx->nominal();
                    c.accrual = fx->accrualPeriod();
                    c.signedAmount = sign * fx->amount();
                    c.sign = sign;
                    bFixed.push_back(c);
                    couponSpans.push_back({c.startIdx, c.payIdx});
                } else if (auto cpn = QuantLib::ext::dynamic_pointer_cast<Coupon>(cf)) {
                    // compounded floats in B stay stochastic: the mapper recorded them; their spans
                    // are checked below through the instrument
                    (void)cpn;
                }
                lastBPay = std::max(lastBPay, grid.times()[grid.index(cf->date(), toleranceDays, what.str() + " flow")]);
            }
        } else if (args.currency[l] == inputs.ccyA) {
            for (const auto& cf : args.legs[l]) {
                if (cf->date() <= today)
                    continue;
                const Size payIdx = grid.index(cf->date(), toleranceDays, what.str() + " pay date");
                lastAPay = std::max(lastAPay, grid.times()[payIdx]);
                auto onQle = QuantLib::ext::dynamic_pointer_cast<QuantExt::OvernightIndexedCoupon>(cf);
                auto onQl = QuantLib::ext::dynamic_pointer_cast<QuantLib::OvernightIndexedCoupon>(cf);
                if (onQle || onQl) {
                    TelescopedCoupon t;
                    if (onQle) {
                        checkOvernightCoupon(*onQle, today, inputs.discountA, what.str());
                        QL_REQUIRE(!onQle->includeSpread(),
                                   "fmmMapXccyReduced: " << what.str() << ": includeSpread compounding is not supported");
                        QL_REQUIRE(onQle->lookback() == 0 * Days && onQle->rateCutoff() == 0,
                                   "fmmMapXccyReduced: " << what.str() << ": lookback / rate cut-off are not supported");
                        t.start = onQle->accrualStartDate();
                        t.nominal = onQle->nominal();
                        t.spread = onQle->spread();
                        t.accrual = onQle->accrualPeriod();
                        t.startIdx = grid.index(t.start, toleranceDays, what.str() + " accrual start");
                        QL_REQUIRE(grid.index(onQle->accrualEndDate(), toleranceDays, what.str() + " accrual end") == payIdx,
                                   "fmmMapXccyReduced: " << what.str() << ": payment lag beyond the grid tolerance");
                    } else {
                        checkOvernightCoupon(*onQl, today, inputs.discountA, what.str());
                        QL_REQUIRE(onQl->lockoutDays() == 0,
                                   "fmmMapXccyReduced: " << what.str() << ": lockout conventions are not supported");
                        t.start = onQl->accrualStartDate();
                        t.nominal = onQl->nominal();
                        t.spread = onQl->spread();
                        t.accrual = onQl->accrualPeriod();
                        t.startIdx = grid.index(t.start, toleranceDays, what.str() + " accrual start");
                        QL_REQUIRE(grid.index(onQl->accrualEndDate(), toleranceDays, what.str() + " accrual end") == payIdx,
                                   "fmmMapXccyReduced: " << what.str() << ": payment lag beyond the grid tolerance");
                    }
                    t.payIdx = payIdx;
                    t.pay = cf->date();
                    t.sign = sign;
                    QL_REQUIRE(t.payIdx > t.startIdx, "fmmMapXccyReduced: " << what.str() << ": inconsistent coupon dates");
                    // rate part: +N at the accrual start (a coupon start: entered when it falls on the
                    // settlement date), -N at the accrual end; spread: N s tau at the payment
                    aStartFlowsB[t.startIdx] += t.sign * t.nominal * fwdFx(t.start);
                    inst.fixedFlows[t.payIdx] -= t.sign * t.nominal * fwdFx(t.pay);
                    const Real spreadAmountB = t.sign * t.nominal * t.spread * t.accrual * fwdFx(t.pay);
                    inst.fixedFlows[t.payIdx] += spreadAmountB;
                    diag.spreadStreamValueB += spreadAmountB * inputs.discountB->discount(t.pay);
                    inst.lastFlowIdx = std::max(inst.lastFlowIdx, t.payIdx);
                    aFloats.push_back(t);
                    couponSpans.push_back({t.startIdx, t.payIdx});
                    ++diag.telescopedCoupons;
                } else if (auto fx = QuantLib::ext::dynamic_pointer_cast<FixedRateCoupon>(cf)) {
                    QL_REQUIRE(fx->accrualStartDate() >= today || fx->accrualEndDate() == fx->date(),
                               "fmmMapXccyReduced: " << what.str() << ": fixed coupon accruing at the valuation date");
                    inst.fixedFlows[payIdx] += sign * fx->amount() * fwdFx(fx->date());
                    inst.lastFlowIdx = std::max(inst.lastFlowIdx, payIdx);
                    const Size s = fx->accrualStartDate() > today
                                       ? grid.index(fx->accrualStartDate(), toleranceDays, what.str() + " accrual start")
                                       : 0;
                    couponSpans.push_back({s, payIdx});
                    ++diag.convertedFlows;
                } else if (QuantLib::ext::dynamic_pointer_cast<Coupon>(cf)) {
                    QL_FAIL("fmmMapXccyReduced: " << what.str() << ": unsupported coupon type in currency A");
                } else {
                    // notional exchanges, redemptions: deterministic A amounts
                    inst.fixedFlows[payIdx] += sign * cf->amount() * fwdFx(cf->date());
                    inst.lastFlowIdx = std::max(inst.lastFlowIdx, payIdx);
                    ++diag.convertedFlows;
                }
            }
        } else {
            QL_FAIL("fmmMapXccyReduced: leg " << l << " currency " << args.currency[l].code() << " is neither "
                                              << inputs.ccyA.code() << " nor " << inputs.ccyB.code());
        }
    }
    for (const auto& f : inst.compoundedFloats)
        couponSpans.push_back({f.startIdx, f.payIdx});
    // coupon starts that do not fall on a settlement date are ordinary flows
    for (Size j = 1; j <= M; ++j) {
        inst.fixedFlows[j] += aStartFlowsB[j];
        if (std::fabs(aStartFlowsB[j]) > 0.0)
            inst.lastFlowIdx = std::max(inst.lastFlowIdx, j);
    }
    QL_REQUIRE(inst.lastFlowIdx >= 1, "fmmMapXccyReduced: no future cash flows");

    usedSettle.clear();
    if (!args.exercise)
        return inst;
    QL_REQUIRE(args.exercise->type() != Exercise::American,
               "fmmMapXccyReduced: American exercise is not supported (grid-date Bermudan only)");
    auto rebated = QuantLib::ext::dynamic_pointer_cast<RebatedExercise>(args.exercise);
    const std::vector<Date>& noticeDates = args.exercise->dates();
    const std::vector<Date> settleDates = args.settlementDates.empty() ? noticeDates : args.settlementDates;
    QL_REQUIRE(settleDates.size() == noticeDates.size(), "fmmMapXccyReduced: " << settleDates.size()
                                                                              << " settlement dates for "
                                                                              << noticeDates.size()
                                                                              << " exercise dates");
    for (Size i = 0; i < noticeDates.size(); ++i) {
        const Date& dn = noticeDates[i];
        const Date& ds = settleDates[i];
        if (dn <= today)
            continue;
        QL_REQUIRE(ds >= dn, "fmmMapXccyReduced: settlement date " << ds << " before notice date " << dn);
        const Size noticeIdx = grid.index(dn, toleranceDays, "exercise (notice) date");
        const Size settleIdx = std::max(noticeIdx, grid.index(ds, toleranceDays, "settlement date"));
        if (settleIdx >= inst.lastFlowIdx)
            continue; // entry at or after the last flow has no value
        if (!inst.rights.empty() && noticeIdx <= inst.rights.back().noticeIdx)
            continue; // duplicates after mapping are collapsed
        for (const auto& span : couponSpans)
            QL_REQUIRE(!(span.first < settleIdx && span.second > settleIdx),
                       "fmmMapXccyReduced: settlement date " << ds << " falls inside a coupon period (grid "
                                                             << grid.dates()[span.first] << " - "
                                                             << grid.dates()[span.second]
                                                             << "); whole-period exercise required");
        FmmXccyReductionRight r;
        r.noticeDate = dn;
        r.settleDate = ds;
        r.noticeIdx = noticeIdx;
        r.settleIdx = settleIdx;
        r.expiry = grid.times()[noticeIdx];
        r.forwardFxAtSettle = fwdFx(ds);
        // the fee: the A coupon start on the settlement date (entered), plus the rebates
        Real fee = aStartFlowsB[settleIdx];
        if (rebated) {
            const Date payDate = rebated->rebatePaymentDate(dn);
            QL_REQUIRE(grid.index(payDate, toleranceDays, "rebate payment date") == settleIdx,
                       "fmmMapXccyReduced: rebate payment date " << payDate << " is not the settlement date " << ds);
            for (Size no = 0; no < rebated->rebateCurrencies().size(); ++no) {
                const Real reb = rebated->rebate(i, no);
                if (reb == 0.0)
                    continue;
                const Currency& c = rebated->rebateCurrency(no);
                if (c == inputs.ccyB)
                    fee += reb;
                else if (c == inputs.ccyA)
                    fee += reb * r.forwardFxAtSettle;
                else
                    QL_FAIL("fmmMapXccyReduced: rebate currency " << c.code() << " is neither A nor B");
            }
        }
        r.feeFlowB = fee;
        // B fixed leg after settlement: notional, annuity, fixed rate, convexity ratio (from the notice date)
        Real annuityB = 0.0, convexityB = 0.0, signedFixed = 0.0, signB = 0.0;
        for (const auto& c : bFixed) {
            if (c.payIdx <= settleIdx)
                continue;
            QL_REQUIRE(r.notionalB == 0.0 || std::fabs(c.nominal - r.notionalB) < 1e-8 * c.nominal,
                       "fmmMapXccyReduced: the B fixed leg must have a constant notional");
            r.notionalB = c.nominal;
            signB = c.sign;
            const Real df = inputs.discountB->discount(c.pay);
            annuityB += c.accrual * df;
            convexityB += c.accrual * df * (grid.times()[c.payIdx] - r.expiry);
            signedFixed += c.signedAmount * df; // annuity-weighted: R_fix for a constant-rate leg
            r.swapLengthB = std::max(r.swapLengthB, grid.times()[c.payIdx] - grid.times()[settleIdx]);
        }
        QL_REQUIRE(annuityB > 0.0 && r.notionalB > 0.0,
                   "fmmMapXccyReduced: no B fixed coupon after the settlement date " << ds);
        r.annuityB = annuityB;
        r.gammaB = convexityB / annuityB;
        r.fixedRate = std::fabs(signedFixed) / (r.notionalB * annuityB);
        // A floating leg after settlement: spread value, annuity, ATM rate, convexity ratio
        Real annuityA = 0.0, convexityA = 0.0, spreadValueB = 0.0;
        Date firstA, lastA;
        for (const auto& t : aFloats) {
            if (t.startIdx < settleIdx)
                continue;
            QL_REQUIRE(r.notionalA == 0.0 || std::fabs(t.nominal - r.notionalA) < 1e-8 * t.nominal,
                       "fmmMapXccyReduced: the A floating leg must have a constant notional");
            r.notionalA = t.nominal;
            const Real dfA = inputs.discountA->discount(t.pay);
            annuityA += t.accrual * dfA;
            convexityA += t.accrual * dfA * (grid.times()[t.payIdx] - r.expiry);
            spreadValueB += t.sign * t.nominal * t.spread * t.accrual * fwdFx(t.pay) * inputs.discountB->discount(t.pay);
            if (firstA == Date() || t.start < firstA)
                firstA = t.start;
            if (lastA == Date() || t.pay > lastA)
                lastA = t.pay;
            r.swapLengthA = std::max(r.swapLengthA, grid.times()[t.payIdx] - grid.times()[settleIdx]);
        }
        if (annuityA > 0.0) {
            r.annuityA = annuityA;
            r.gammaA = convexityA / annuityA;
            r.forwardSwapRateA = (inputs.discountA->discount(firstA) - inputs.discountA->discount(lastA)) / annuityA;
        }
        // strike (eq. 3.3a): the fixed rate net of the equivalent spread, i.e. the fixed leg's rate
        // moved by the spread stream's B value per unit of B annuity, on the fixed leg's side
        const Real sEqHolder = spreadValueB / (r.notionalB * annuityB);
        r.strike = r.fixedRate + signB * sEqHolder;
        r.equivalentSpread = r.fixedRate - r.strike;
        inst.rights.push_back({noticeIdx, settleIdx, fee});
        usedSettle.push_back(ds);
        diag.rights.push_back(r);
    }
    (void)lastAPay;
    (void)lastBPay;
    inst.validate(M, false);
    return inst;
}

} // namespace QuantExt
