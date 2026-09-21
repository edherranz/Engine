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
#include <qle/pricingengines/fmmlsmmultilegoptionengine.hpp>

#include <ql/cashflows/capflooredcoupon.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/cashflows/iborcoupon.hpp>
#include <ql/cashflows/overnightindexedcoupon.hpp>
#include <ql/cashflows/simplecashflow.hpp>

#include <cmath>
#include <memory>
#include <sstream>

namespace QuantExt {

FmmLegMapper::FmmLegMapper(const FmmGrid& grid, const Natural toleranceDays, const Date& today)
    : grid_(grid), toleranceDays_(toleranceDays), today_(today) {}

void FmmLegMapper::addLeg(const Leg& leg, const Real sign, FmmCallableInstrument& inst, const std::string& what) const {
    const Size M = grid_.numberOfRates();
    if (inst.fixedFlows.size() != M + 1) {
        inst.fixedFlows.assign(M + 1, 0.0);
        inst.floatWeights.assign(M + 1, 0.0);
    }
    for (const auto& cf : leg) {
        if (cf->date() <= today_)
            continue; // already paid (or paying today): not part of the forward value
        if (auto c = QuantLib::ext::dynamic_pointer_cast<CappedFlooredCoupon>(cf)) {
            QL_FAIL("FmmLegMapper: " << what << ": capped/floored coupons are not supported in this release");
        } else if (auto on = QuantLib::ext::dynamic_pointer_cast<QuantExt::OvernightIndexedCoupon>(cf)) {
            QL_REQUIRE(on->accrualStartDate() >= today_,
                       "FmmLegMapper: " << what << ": coupon accruing at the valuation date (" << on->accrualStartDate()
                                        << " - " << on->accrualEndDate() << ") is not supported (already-started underlying)");
            QL_REQUIRE(!on->includeSpread(), "FmmLegMapper: " << what << ": includeSpread compounding is not supported");
            QL_REQUIRE(on->lookback() == 0 * Days && on->rateCutoff() == 0,
                       "FmmLegMapper: " << what
                                        << ": lookback / rate cut-off conventions are not supported in this "
                                           "release (baseline scope: plain compounded RFR)");
            // the rate window must be the accrual period (up to calendar adjustment): ORE's fixing-in-advance
            // (IsInArrears=false) and last-recent-period legs shift the value dates by whole periods
            QL_REQUIRE(!on->valueDates().empty(), "FmmLegMapper: " << what << ": overnight coupon without value dates");
            const auto dStart = on->valueDates().front() - on->accrualStartDate();
            const auto dEnd = on->valueDates().back() - on->accrualEndDate();
            QL_REQUIRE(dStart >= -7 && dStart <= 7 && dEnd >= -7 && dEnd <= 7,
                       "FmmLegMapper: " << what << ": the rate computation window (" << on->valueDates().front()
                                        << " - " << on->valueDates().back() << ") differs from the accrual period ("
                                        << on->accrualStartDate() << " - " << on->accrualEndDate()
                                        << "): fixing in advance / last recent period is not supported (baseline "
                                           "scope: in-arrears compounded RFR)");
            FmmCallableInstrument::CompoundedFloat f;
            f.startIdx = grid_.index(on->accrualStartDate(), toleranceDays_, what + " accrual start");
            f.endIdx = grid_.index(on->accrualEndDate(), toleranceDays_, what + " accrual end");
            f.payIdx = grid_.index(on->date(), toleranceDays_, what + " pay date");
            QL_REQUIRE(f.endIdx > f.startIdx && f.payIdx >= f.endIdx, "FmmLegMapper: " << what << ": inconsistent coupon dates");
            f.weight = sign * on->nominal() * on->gearing();
            f.spreadAmount = sign * on->nominal() * on->spread() * on->accrualPeriod();
            inst.compoundedFloats.push_back(f);
            inst.lastFlowIdx = std::max(inst.lastFlowIdx, f.payIdx);
        } else if (auto onq = QuantLib::ext::dynamic_pointer_cast<QuantLib::OvernightIndexedCoupon>(cf)) {
            QL_REQUIRE(onq->accrualStartDate() >= today_,
                       "FmmLegMapper: " << what << ": coupon accruing at the valuation date is not supported");
            QL_REQUIRE(onq->lockoutDays() == 0,
                       "FmmLegMapper: " << what << ": lockout conventions are not supported in this release");
            FmmCallableInstrument::CompoundedFloat f;
            f.startIdx = grid_.index(onq->accrualStartDate(), toleranceDays_, what + " accrual start");
            f.endIdx = grid_.index(onq->accrualEndDate(), toleranceDays_, what + " accrual end");
            f.payIdx = grid_.index(onq->date(), toleranceDays_, what + " pay date");
            QL_REQUIRE(f.endIdx > f.startIdx && f.payIdx >= f.endIdx, "FmmLegMapper: " << what << ": inconsistent coupon dates");
            f.weight = sign * onq->nominal() * onq->gearing();
            f.spreadAmount = sign * onq->nominal() * onq->spread() * onq->accrualPeriod();
            inst.compoundedFloats.push_back(f);
            inst.lastFlowIdx = std::max(inst.lastFlowIdx, f.payIdx);
        } else if (auto ib = QuantLib::ext::dynamic_pointer_cast<IborCoupon>(cf)) {
            QL_FAIL("FmmLegMapper: " << what
                                     << ": term-rate (ibor) coupons are outside the single-currency RFR baseline of "
                                        "this release");
        } else if (auto fx = QuantLib::ext::dynamic_pointer_cast<FixedRateCoupon>(cf)) {
            QL_REQUIRE(fx->accrualStartDate() >= today_ || fx->accrualEndDate() == fx->date(),
                       "FmmLegMapper: " << what << ": fixed coupon accruing at the valuation date");
            const Size idx = grid_.index(fx->date(), toleranceDays_, what + " fixed pay date");
            inst.fixedFlows[idx] += sign * fx->amount();
            inst.lastFlowIdx = std::max(inst.lastFlowIdx, idx);
        } else if (QuantLib::ext::dynamic_pointer_cast<Coupon>(cf)) {
            QL_FAIL("FmmLegMapper: " << what << ": unsupported coupon type");
        } else {
            // notional exchanges, redemptions, fees
            const Size idx = grid_.index(cf->date(), toleranceDays_, what + " cash flow date");
            inst.fixedFlows[idx] += sign * cf->amount();
            inst.lastFlowIdx = std::max(inst.lastFlowIdx, idx);
        }
    }
}

Real FmmLegMapper::underlyingCurveValue(const FmmCallableInstrument& inst, const FmmParametrization& p) {
    return fmmUnderlyingCurveValue(inst, p);
}

void fmmWriteLsmResults(std::map<std::string, QuantLib::ext::any>& out, const FmmLsmResult& res,
                        const FmmDualBoundResult* dual, const FmmLsmEngineConfig& config, const FmmGrid& grid) {
    out["fmmLsmLowerBound"] = res.lowerBound;
    out["fmmLsmLowerBoundStdError"] = res.lowerBoundSe;
    out["fmmLsmTrainingValue"] = res.trainingValue;
    out["fmmLsmTrainingPaths"] = static_cast<Real>(config.lsm.trainingPaths);
    out["fmmLsmValuationPaths"] = static_cast<Real>(config.lsm.valuationPaths);
    out["fmmLsmTrainingSeed"] = static_cast<Real>(config.lsm.trainingSeed);
    out["fmmLsmValuationSeed"] = static_cast<Real>(config.lsm.valuationSeed);
    out["fmmLsmBasisOrder"] = static_cast<Real>(config.lsm.basisOrder);
    out["fmmNoExerciseProbability"] = res.noExerciseProbability;
    out["fmmExpectedExerciseTime"] = res.expectedExerciseTime;
    out["fmmExerciseProbabilityType"] = std::string("unconditional");
    std::vector<Real> probs = res.exerciseProbability, probSe = res.exerciseProbabilitySe;
    out["fmmExerciseProbabilities"] = probs;
    out["fmmExerciseProbabilityStdErrors"] = probSe;
    if (dual) {
        out["fmmDualityGap"] = dual->gap;
        out["fmmDualityGapStdError"] = dual->gapSe;
        out["fmmLsmUpperBound"] = res.lowerBound + dual->gap;
        out["fmmLsmUpperBoundStdError"] = std::sqrt(res.lowerBoundSe * res.lowerBoundSe + dual->gapSe * dual->gapSe);
        out["fmmDualOuterPaths"] = static_cast<Real>(dual->outerPaths);
        out["fmmDualInnerPaths"] = static_cast<Real>(dual->innerPaths);
    }
    out["fmmGridPeriods"] = static_cast<Real>(grid.numberOfRates());
    out["fmmGridMaxDateMismatchDays"] = static_cast<Real>(grid.maxMismatchDays());
    out["fmmRuntimeSeconds"] = res.runtimeSeconds + (dual ? dual->runtimeSeconds : 0.0);
}

FmmLsmMultiLegOptionEngine::FmmLsmMultiLegOptionEngine(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                                                       const QuantLib::ext::shared_ptr<FmmGrid>& grid,
                                                       const FmmLsmEngineConfig& config)
    : model_(model), grid_(grid), config_(config) {
    QL_REQUIRE(model_ && grid_, "FmmLsmMultiLegOptionEngine: model or grid is null");
    QL_REQUIRE(grid_->numberOfRates() == model_->parametrization()->numberOfRates(),
               "FmmLsmMultiLegOptionEngine: grid and model disagree on the number of periods");
    registerWith(model_->parametrization()->termStructure());
    for (const auto& o : config_.observables)
        registerWith(o);
}

FmmCallableInstrument fmmMapLegs(const std::vector<Leg>& legs, const std::vector<Real>& signs, const FmmGrid& grid,
                                 const Natural toleranceDays, const Date& today) {
    QL_REQUIRE(legs.size() == signs.size(), "fmmMapLegs: " << legs.size() << " legs, " << signs.size() << " signs");
    FmmCallableInstrument inst;
    inst.style = FmmCallableInstrument::Style::Cancel;
    FmmLegMapper mapper(grid, toleranceDays, today);
    for (Size i = 0; i < legs.size(); ++i) {
        std::ostringstream what;
        what << "leg " << i;
        mapper.addLeg(legs[i], signs[i] < 0.0 ? -1.0 : 1.0, inst, what.str());
    }
    QL_REQUIRE(inst.lastFlowIdx >= 1, "fmmMapLegs: no future cash flows");
    return inst;
}

FmmCallableInstrument fmmMapMultiLegOption(const MultiLegOption::arguments& args, const FmmGrid& grid,
                                           const Natural toleranceDays, const Date& today,
                                           std::vector<Date>& usedSettle) {
    QL_REQUIRE(args.legs.size() == args.payer.size(), "fmmMapMultiLegOption: inconsistent legs / payer vectors");
    QL_REQUIRE(args.settlementType == Settlement::Physical,
               "fmmMapMultiLegOption: physical settlement only (cash settlement is not supported)");
    QL_REQUIRE(!args.midCouponExercise, "fmmMapMultiLegOption: mid-coupon exercise is not supported");
    // flows from the OPTION HOLDER's view: receive legs with payer = false, pay legs with payer = true
    std::vector<Real> signs;
    for (const bool payer : args.payer)
        signs.push_back(payer ? -1.0 : 1.0);
    FmmCallableInstrument inst = fmmMapLegs(args.legs, signs, grid, toleranceDays, today);
    inst.style = FmmCallableInstrument::Style::Enter;
    usedSettle.clear();
    if (!args.exercise)
        return inst; // no exercise: the instrument is the underlying itself
    QL_REQUIRE(args.exercise->type() != Exercise::American,
               "fmmMapMultiLegOption: American exercise is not supported (grid-date Bermudan only)");
    // exercise dates -> rights; every exercise date must be a grid date at which all legs have an
    // accrual boundary, i.e. no coupon straddles it (checked through the flow mapping: a coupon
    // with accrual start before and pay date after an exercise date would be mis-assigned)
    // ORE passes the notice dates as the exercise dates and the settlement (swap entry) dates
    // alongside: a right decides at the notice grid date and enters the flows after the settlement
    // grid date (A4 notice-period rights, adapted payoff in the dual bound)
    const std::vector<Date>& noticeDates = args.exercise->dates();
    const std::vector<Date> settleDates = args.settlementDates.empty() ? noticeDates : args.settlementDates;
    QL_REQUIRE(settleDates.size() == noticeDates.size(), "fmmMapMultiLegOption: " << settleDates.size()
                                                                                  << " settlement dates for "
                                                                                  << noticeDates.size()
                                                                                  << " exercise dates");
    for (Size i = 0; i < noticeDates.size(); ++i) {
        const Date& dn = noticeDates[i];
        const Date& ds = settleDates[i];
        if (dn <= today)
            continue;
        QL_REQUIRE(ds >= dn, "fmmMapMultiLegOption: settlement date " << ds << " before notice date " << dn);
        const Size noticeIdx = grid.index(dn, toleranceDays, "exercise (notice) date");
        const Size settleIdx = std::max(noticeIdx, grid.index(ds, toleranceDays, "settlement date"));
        if (settleIdx >= inst.lastFlowIdx)
            continue; // entry at or after the last flow has no value
        for (const auto& f : inst.compoundedFloats)
            QL_REQUIRE(!(f.startIdx < settleIdx && f.payIdx > settleIdx),
                       "fmmMapMultiLegOption: settlement date "
                           << ds << " falls inside a compounded coupon period (" << grid.dates()[f.startIdx] << " - "
                           << grid.dates()[f.endIdx] << "); whole-period exercise required");
        inst.rights.push_back({noticeIdx, settleIdx, 0.0});
        usedSettle.push_back(ds);
    }
    // notice indices strictly ascending (duplicates after mapping are collapsed)
    std::vector<FmmCallableInstrument::Right> rights;
    for (const auto& r : inst.rights)
        if (rights.empty() || r.noticeIdx > rights.back().noticeIdx)
            rights.push_back(r);
    inst.rights = rights;
    return inst;
}

void FmmLsmMultiLegOptionEngine::calculate() const {
    const auto& p = *model_->parametrization();
    const Date today = grid_->referenceDate();
    QL_REQUIRE(p.termStructure()->referenceDate() == today,
               "FmmLsmMultiLegOptionEngine: model curve reference date " << p.termStructure()->referenceDate()
                                                                         << " != grid reference date " << today);
    QL_REQUIRE(arguments_.legs.size() == arguments_.payer.size() && arguments_.legs.size() == arguments_.currency.size(),
               "FmmLsmMultiLegOptionEngine: inconsistent legs / payer / currency vectors");
    for (const auto& c : arguments_.currency)
        QL_REQUIRE(c == p.currency(), "FmmLsmMultiLegOptionEngine: leg currency " << c << " != model currency "
                                                                                  << p.currency());

    std::vector<Date> usedSettle;
    FmmCallableInstrument inst = fmmMapMultiLegOption(arguments_, *grid_, config_.gridToleranceDays, today, usedSettle);
    const Real underlyingNpv = FmmLegMapper::underlyingCurveValue(inst, p);
    results_.underlyingNpv = underlyingNpv;
    results_.additionalResults["underlyingNpv"] = underlyingNpv;
    if (inst.rights.empty()) {
        // no exercise object: the instrument is the underlying itself; no right left: worthless
        results_.value = arguments_.exercise ? 0.0 : underlyingNpv;
        results_.additionalResults["fmmGridMaxDateMismatchDays"] = static_cast<Real>(grid_->maxMismatchDays());
        return;
    }

    FmmLsmPricer pricer(model_, inst, config_.lsm);
    FmmLsmResult res;
    bool reusedPolicy = false;
    if (config_.policyMode == FmmPolicyMode::Frozen && frozenPolicy_) {
        // frozen-policy risk: the policy trained by the first valuation, revalued on the same paths
        res = pricer.valueWithPolicy(*frozenPolicy_, config_.lsm.valuationSeed);
        reusedPolicy = true;
    } else {
        res = pricer.calculate();
        if (config_.policyMode == FmmPolicyMode::Frozen)
            frozenPolicy_ = pricer.policy();
    }
    results_.value = res.lowerBound;
    std::unique_ptr<FmmDualBoundResult> dual;
    if (config_.dualBound && !reusedPolicy)
        dual = std::make_unique<FmmDualBoundResult>(
            pricer.dualBound(config_.dualOuterPaths, config_.dualInnerPaths, config_.dualSeed));
    fmmWriteLsmResults(results_.additionalResults, res, dual.get(), config_, *grid_);
    results_.additionalResults["fmmExerciseDates"] = arguments_.exercise->dates();
    results_.additionalResults["fmmPolicyMode"] =
        std::string(config_.policyMode == FmmPolicyMode::Frozen
                        ? (reusedPolicy ? "Frozen (policy of the first valuation reused)" : "Frozen (policy trained)")
                        : "Retrain");
    results_.additionalResults["fmmSettlementDates"] = usedSettle;
    if (config_.calibrationResults)
        for (const auto& kv : config_.calibrationResults())
            results_.additionalResults[kv.first] = kv.second;
}

} // namespace QuantExt
