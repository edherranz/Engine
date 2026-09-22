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

#include "toplevelfixture.hpp"
#include <boost/test/unit_test.hpp>

#include <qle/instruments/multilegoption.hpp>
#include <qle/instruments/rebatedexercise.hpp>
#include <qle/models/fmmanalytics.hpp>
#include <qle/models/fmmgrid.hpp>
#include <qle/models/fmmlsmpricer.hpp>
#include <qle/models/fmmxccyreduction.hpp>
#include <qle/models/forwardmarketmodel.hpp>
#include <qle/pricingengines/fmmlsmmultilegoptionengine.hpp>
#include <qle/pricingengines/fmmlsmxccyreducedengine.hpp>

#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/cashflows/overnightindexedcoupon.hpp>
#include <ql/cashflows/simplecashflow.hpp>
#include <ql/currencies/america.hpp>
#include <ql/currencies/europe.hpp>
#include <ql/exercise.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/schedule.hpp>

#include <cmath>
#include <set>

using namespace QuantLib;
using namespace QuantExt;

BOOST_FIXTURE_TEST_SUITE(QuantExtTestSuite, qle::test::TopLevelFixture)

BOOST_AUTO_TEST_SUITE(FmmXccyTest)

namespace {

/* A9.1 test-bed: a cancellable cross-currency swap from the holder's view of the cancellation
   right, i.e. the Bermudan option into the reverse swap (Hagan eqs. 2.4-2.5): receive USD overnight
   plus a spread on N_A with the final exchange, pay EUR fixed on N_B with the final exchange,
   rebates on the exercise date: pay N_A (USD), receive N_B (EUR). USD flat 4%, EUR flat 2.5%,
   EUR/USD spot 1.03 USD per EUR (B = EUR per unit of A = USD: 1/1.03), quarterly USD and annual
   EUR schedules starting three months after the valuation date, five years, cancellable yearly. */
struct XccyBed {
    XccyBed(const Real spread, const bool rebatesEqualNotionals = true)
        : asof(19, September, 2026), dc(Actual365Fixed()), spread(spread) {
        Settings::instance().evaluationDate() = asof;
        curveA = Handle<YieldTermStructure>(QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.04, dc));
        curveB = Handle<YieldTermStructure>(QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.025, dc));
        fxBperA = Handle<Quote>(QuantLib::ext::make_shared<SimpleQuote>(1.0 / 1.03));
        start = asof + Period(3, Months);
        end = start + Period(5, Years);
        schedA = Schedule(start, end, Period(3, Months), NullCalendar(), Unadjusted, Unadjusted, DateGeneration::Forward, false);
        schedB = Schedule(start, end, Period(1, Years), NullCalendar(), Unadjusted, Unadjusted, DateGeneration::Forward, false);
        indexA = QuantLib::ext::make_shared<OvernightIndex>("USDON", 0, USDCurrency(), NullCalendar(), dc, curveA);
        indexB = QuantLib::ext::make_shared<OvernightIndex>("EURON", 0, EURCurrency(), NullCalendar(), dc, curveB);
        legA = OvernightLeg(schedA, indexA).withNotionals(NA).withSpreads(spread);
        legA.push_back(QuantLib::ext::make_shared<SimpleCashFlow>(NA, end));
        legB = FixedRateLeg(schedB).withNotionals(NB).withCouponRates(R, dc);
        legB.push_back(QuantLib::ext::make_shared<SimpleCashFlow>(NB, end));
        for (Size k = 1; k < schedB.size() - 1; ++k)
            exDates.push_back(schedB[k]);
        const Real rebA = rebatesEqualNotionals ? -NA : -0.5 * NA;
        std::vector<std::vector<Real>> rebates(exDates.size(), std::vector<Real>{rebA, NB});
        exercise = QuantLib::ext::make_shared<RebatedExercise>(BermudanExercise(exDates), exDates, rebates,
                                                               std::vector<Currency>{USDCurrency(), EURCurrency()},
                                                               0 * Days);
        std::set<Date> dates(schedA.dates().begin(), schedA.dates().end());
        dates.insert(schedB.dates().begin(), schedB.dates().end());
        grid = QuantLib::ext::make_shared<FmmGrid>(asof, dates, Period(3, Months));
        const Size M = grid->numberOfRates();
        Array shifts(M);
        for (Size j = 1; j <= M; ++j)
            shifts[j - 1] = 1.0 / (grid->times()[j] - grid->times()[j - 1]);
        std::vector<Array> volLevels(M, Array(1, 0.0025));
        param = QuantLib::ext::make_shared<FmmParametrization>(EURCurrency(), curveB, grid->times(), shifts, Array(),
                                                               volLevels,
                                                               FmmParametrization::LocalVolType::DisplacedDiffusion, 0.6,
                                                               0.08, 3);
        model = QuantLib::ext::make_shared<ForwardMarketModel>(param);
        inputs.ccyA = USDCurrency();
        inputs.ccyB = EURCurrency();
        inputs.discountA = curveA;
        inputs.discountB = curveB;
        inputs.fxSpotBperA = fxBperA;
    }
    MultiLegOption::arguments arguments(const bool withExercise) const {
        auto option = QuantLib::ext::make_shared<MultiLegOption>(
            std::vector<Leg>{legA, legB}, std::vector<bool>{false, true},
            std::vector<Currency>{USDCurrency(), EURCurrency()},
            withExercise ? exercise : QuantLib::ext::shared_ptr<Exercise>(), Settlement::Physical,
            Settlement::PhysicalOTC, withExercise ? exDates : std::vector<Date>());
        MultiLegOption::arguments args;
        option->setupArguments(&args);
        return args;
    }
    Real fwdFx(const Date& d) const { return fxBperA->value() * curveA->discount(d) / curveB->discount(d); }

    Date asof, start, end;
    Actual365Fixed dc;
    Real spread;
    const Real NA = 10000000.0, NB = 9700000.0, R = 0.0225;
    Handle<YieldTermStructure> curveA, curveB;
    Handle<Quote> fxBperA;
    Schedule schedA, schedB;
    QuantLib::ext::shared_ptr<OvernightIndex> indexA, indexB;
    Leg legA, legB;
    std::vector<Date> exDates;
    QuantLib::ext::shared_ptr<Exercise> exercise;
    QuantLib::ext::shared_ptr<FmmGrid> grid;
    QuantLib::ext::shared_ptr<FmmParametrization> param;
    QuantLib::ext::shared_ptr<ForwardMarketModel> model;
    FmmXccyReductionInputs inputs;
};

} // namespace

BOOST_AUTO_TEST_CASE(testSwapSpecBasisFlows) {
    BOOST_TEST_MESSAGE("A9.1: deterministic basis flows on the floating leg of an FmmSwapSpec: the forward swap rate "
                       "is (P_a - P_b + sum c_j P_j) / A, the approximation's ATM implied vol agrees with the "
                       "same-model Monte Carlo, and zero basis reproduces the plain spec...");
    XccyBed bed(0.0);
    const auto& p = *bed.param;
    FmmSwapSpec plain;
    plain.a = 4;
    plain.b = 20;
    for (Size j = 8; j <= 20; j += 4) {
        plain.fixedPayIndices.push_back(j);
        plain.fixedAccruals.push_back(p.rateTime(j) - p.rateTime(j - 4));
    }
    FmmSwapSpec withBasis = plain;
    // a 15 bp deterministic basis paid quarterly on the floating leg (times the accrual)
    Real basisPv = 0.0;
    for (Size j = 5; j <= 20; ++j) {
        withBasis.basisPayIndices.push_back(j);
        withBasis.basisAmounts.push_back(0.0015 * p.tau(j));
        basisPv += 0.0015 * p.tau(j) * bed.curveB->discount(p.rateTime(j));
    }
    const Real A = fmmAnnuity(p, plain);
    const Real S0 = fmmForwardSwapRate(p, plain);
    const Real S1 = fmmForwardSwapRate(p, withBasis);
    BOOST_TEST_MESSAGE("forward swap rate " << S0 << " -> with basis " << S1 << " (basis PV / annuity " << basisPv / A << ")");
    BOOST_CHECK_SMALL(S1 - (S0 + basisPv / A), 1e-14);
    BOOST_CHECK_SMALL(fmmAnnuity(p, withBasis) - A, 1e-14);
    for (const auto method : {FmmSwaptionApproxMethod::StrikeFlatNormal, FmmSwaptionApproxMethod::EffectiveShift}) {
        const auto r0 = fmmSwaptionApprox(p, plain, S0, Option::Call, method);
        const auto r1 = fmmSwaptionApprox(p, withBasis, S1, Option::Call, method);
        const auto mc = fmmSwaptionMc(bed.model, {withBasis}, {S1}, Option::Call, 16384, 8, 4711).front();
        BOOST_TEST_MESSAGE("ATM normal vol: plain " << r0.normalVol * 1e4 << " bp, with basis " << r1.normalVol * 1e4
                                                     << " bp, MC with basis " << mc.normalVol * 1e4 << " +/- "
                                                     << mc.normalVolSe * 1e4 << " bp");
        // the basis moves the forward, not the variance of the rates: the vol changes by well under a bp
        BOOST_CHECK_SMALL(r1.normalVol - r0.normalVol, 1e-4);
        BOOST_CHECK_SMALL(r1.normalVol - mc.normalVol, 1.5e-4 + 3.0 * mc.normalVolSe);
    }
}

BOOST_AUTO_TEST_CASE(testReductionFlowsAndSpread) {
    BOOST_TEST_MESSAGE("A9.1: Hagan's reduction of the cancellation right: with zero spread and rebates equal to the "
                       "notionals the USD side telescopes away exactly (the reduced instrument is the EUR structure "
                       "with the EUR notional as the fee); a spread adds the forward-FX-converted spread stream; the "
                       "equivalent spread and strike follow eqs. 3.3a-b; the reduced underlying equals the two-curve "
                       "value of the swap...");
    // (a) zero spread: only the EUR leg and the EUR rebate remain
    XccyBed bed0(0.0);
    FmmXccyReductionDiagnostics diag0;
    std::vector<Date> settle0;
    const auto args0 = bed0.arguments(true);
    const FmmCallableInstrument inst0 = fmmMapXccyReduced(args0, *bed0.grid, 3, bed0.asof, bed0.inputs, diag0, settle0);
    FmmCallableInstrument eurOnly = fmmMapLegs({bed0.legB}, {-1.0}, *bed0.grid, 3, bed0.asof);
    BOOST_REQUIRE_EQUAL(inst0.fixedFlows.size(), eurOnly.fixedFlows.size());
    // the USD leg telescopes to its forward-start value N_A at the first accrual start (a coupon
    // start before every exercise date, never entered by a right); everything else is the EUR leg
    const Size startIdx = bed0.grid->index(bed0.start, 0, "start");
    for (Size j = 0; j < inst0.fixedFlows.size(); ++j) {
        const Real usdStart = j == startIdx ? bed0.NA * bed0.fwdFx(bed0.start) : 0.0;
        BOOST_CHECK_SMALL(inst0.fixedFlows[j] - eurOnly.fixedFlows[j] - usdStart, 1e-6);
    }
    BOOST_CHECK(inst0.compoundedFloats.empty());
    BOOST_REQUIRE_EQUAL(inst0.rights.size(), bed0.exDates.size());
    BOOST_REQUIRE_EQUAL(diag0.rights.size(), bed0.exDates.size());
    BOOST_CHECK_EQUAL(diag0.telescopedCoupons, bed0.schedA.size() - 1);
    for (Size k = 0; k < inst0.rights.size(); ++k) {
        BOOST_CHECK_SMALL(inst0.rights[k].feeFlow - bed0.NB, 1e-6);
        BOOST_CHECK_SMALL(diag0.rights[k].equivalentSpread, 1e-14);
        BOOST_CHECK_SMALL(diag0.rights[k].strike - bed0.R, 1e-12);
        BOOST_CHECK_SMALL(diag0.rights[k].fixedRate - bed0.R, 1e-12);
        BOOST_CHECK_EQUAL(inst0.rights[k].noticeIdx, inst0.rights[k].settleIdx);
        BOOST_CHECK_EQUAL(bed0.grid->dates()[inst0.rights[k].settleIdx], bed0.exDates[k]);
    }
    // (b) a 20 bp spread: the converted spread stream is the only difference; s_eq per eq. 3.3b
    XccyBed bed(0.0020);
    FmmXccyReductionDiagnostics diag;
    std::vector<Date> settle;
    const auto args = bed.arguments(true);
    const FmmCallableInstrument inst = fmmMapXccyReduced(args, *bed.grid, 3, bed.asof, bed.inputs, diag, settle);
    std::vector<Real> expectedSpread(inst.fixedFlows.size(), 0.0);
    for (Size i = 1; i < bed.schedA.size(); ++i) {
        const Real accrual = bed.dc.yearFraction(bed.schedA[i - 1], bed.schedA[i]);
        expectedSpread[bed.grid->index(bed.schedA[i], 0, "pay")] += bed.NA * 0.0020 * accrual * bed.fwdFx(bed.schedA[i]);
    }
    for (Size j = 0; j < inst.fixedFlows.size(); ++j)
        BOOST_CHECK_SMALL(inst.fixedFlows[j] - inst0.fixedFlows[j] - expectedSpread[j], 1e-6);
    for (Size k = 0; k < diag.rights.size(); ++k) {
        const Size settleIdx = diag.rights[k].settleIdx;
        Real spreadValueB = 0.0, annuityB = 0.0;
        for (Size i = 1; i < bed.schedA.size(); ++i)
            if (bed.grid->index(bed.schedA[i - 1], 0, "start") >= settleIdx)
                spreadValueB += bed.NA * 0.0020 * bed.dc.yearFraction(bed.schedA[i - 1], bed.schedA[i]) *
                                bed.fwdFx(bed.schedA[i]) * bed.curveB->discount(bed.schedA[i]);
        for (Size i = 1; i < bed.schedB.size(); ++i)
            if (bed.grid->index(bed.schedB[i], 0, "pay") > settleIdx)
                annuityB += bed.dc.yearFraction(bed.schedB[i - 1], bed.schedB[i]) * bed.curveB->discount(bed.schedB[i]);
        const Real sEq = spreadValueB / (bed.NB * annuityB);
        BOOST_TEST_MESSAGE("right " << k << " (" << diag.rights[k].noticeDate << "): s_eq " << sEq * 1e4 << " bp, strike "
                                    << diag.rights[k].strike << ", gammaA " << diag.rights[k].gammaA << ", gammaB "
                                    << diag.rights[k].gammaB << ", fee " << diag.rights[k].feeFlowB);
        BOOST_CHECK_SMALL(diag.rights[k].equivalentSpread - sEq, 1e-12);
        BOOST_CHECK_SMALL(diag.rights[k].strike - (bed.R - sEq), 1e-12);
        BOOST_CHECK_SMALL(diag.rights[k].annuityB - annuityB, 1e-12);
        BOOST_CHECK(diag.rights[k].gammaA > 0.0 && diag.rights[k].gammaB > 0.0);
        BOOST_CHECK(diag.rights[k].gammaA <= diag.rights[k].swapLengthA + 1e-12 &&
                    diag.rights[k].gammaB <= diag.rights[k].swapLengthB + 1e-12);
    }
    // (c) the reduced underlying (no exercise) is the swap's two-curve value in EUR: the USD leg on
    // its own curve (par float plus spread plus notional) converted at spot
    const auto argsU = bed.arguments(false);
    FmmXccyReductionDiagnostics diagU;
    std::vector<Date> settleU;
    const FmmCallableInstrument instU = fmmMapXccyReduced(argsU, *bed.grid, 3, bed.asof, bed.inputs, diagU, settleU);
    BOOST_CHECK(instU.rights.empty());
    Real usdLeg = bed.NA * bed.curveA->discount(bed.start); // par float + notional from the (future) start
    for (Size i = 1; i < bed.schedA.size(); ++i)
        usdLeg += bed.NA * 0.0020 * bed.dc.yearFraction(bed.schedA[i - 1], bed.schedA[i]) * bed.curveA->discount(bed.schedA[i]);
    Real eurLeg = bed.NB * bed.curveB->discount(bed.end);
    for (Size i = 1; i < bed.schedB.size(); ++i)
        eurLeg += bed.NB * bed.R * bed.dc.yearFraction(bed.schedB[i - 1], bed.schedB[i]) * bed.curveB->discount(bed.schedB[i]);
    const Real direct = usdLeg * bed.fxBperA->value() - eurLeg;
    const Real reduced = fmmUnderlyingCurveValue(instU, *bed.param);
    BOOST_TEST_MESSAGE("reduced underlying " << reduced << " EUR vs two-curve swap value " << direct << " EUR");
    BOOST_CHECK_SMALL(reduced - direct, 1e-6 * bed.NB);
    // (d) rebates below the notional leave the residual in the fee (converted at the settlement forward)
    XccyBed bedR(0.0020, false);
    FmmXccyReductionDiagnostics diagR;
    std::vector<Date> settleR;
    const FmmCallableInstrument instR = fmmMapXccyReduced(bedR.arguments(true), *bedR.grid, 3, bedR.asof, bedR.inputs, diagR, settleR);
    for (Size k = 0; k < instR.rights.size(); ++k)
        BOOST_CHECK_SMALL(instR.rights[k].feeFlow - (bedR.NB + 0.5 * bedR.NA * bedR.fwdFx(bedR.exDates[k])), 1e-6);
}

BOOST_AUTO_TEST_CASE(testReducedEngineZeroSpreadLimit) {
    BOOST_TEST_MESSAGE("A9.1: with zero spread the reduced engine reproduces the single-currency FMM engine on the "
                       "equivalent EUR Bermudan payer swaption (fixed against the EUR overnight leg) on the same "
                       "grid, model and seeds (Hagan section 3.1; the bank-account identity makes the EUR notional "
                       "flows and the compounded EUR floats pathwise identical)...");
    XccyBed bed(0.0);
    FmmLsmEngineConfig cfg;
    cfg.lsm.trainingPaths = 8192;
    cfg.lsm.valuationPaths = 8192;
    cfg.gridToleranceDays = 3;
    auto option = QuantLib::ext::make_shared<MultiLegOption>(
        std::vector<Leg>{bed.legA, bed.legB}, std::vector<bool>{false, true},
        std::vector<Currency>{USDCurrency(), EURCurrency()}, bed.exercise, Settlement::Physical, Settlement::PhysicalOTC,
        bed.exDates);
    Size gridCalls = 0, modelCalls = 0;
    auto gridProvider = [&](const std::set<Date>& dates) {
        ++gridCalls;
        for (const Date& d : dates)
            BOOST_CHECK_NO_THROW(bed.grid->index(d, 0, "contractual date"));
        return bed.grid;
    };
    auto modelProvider = [&](const FmmXccyReductionDiagnostics& diag, const FmmGrid& grid) {
        ++modelCalls;
        BOOST_CHECK_EQUAL(grid.numberOfRates(), bed.grid->numberOfRates());
        FmmXccyCalibrationRecord rec;
        rec.method = "Standard";
        rec.model = bed.model;
        for (const auto& r : diag.rights) {
            rec.strikes.push_back(r.strike);
            rec.marketVolsB.push_back(0.0);
            rec.targetVols.push_back(0.0);
        }
        return rec;
    };
    option->setPricingEngine(QuantLib::ext::make_shared<FmmLsmXccyReducedEngine>(bed.inputs, EURCurrency(), gridProvider,
                                                                                  modelProvider, cfg));
    const Real reduced = option->NPV();
    BOOST_CHECK_EQUAL(gridCalls, Size(1));
    BOOST_CHECK_EQUAL(modelCalls, Size(1));
    BOOST_CHECK_EQUAL(option->result<std::string>("fmmXccyMethod"), "Standard");
    const auto strikes = option->result<std::vector<Real>>("fmmXccyStrikes");
    BOOST_REQUIRE_EQUAL(strikes.size(), bed.exDates.size());
    for (const Real K : strikes)
        BOOST_CHECK_SMALL(K - bed.R, 1e-12);
    // the single-currency equivalent: pay EUR fixed, receive EUR overnight (quarterly, the USD dates)
    Leg floatB = OvernightLeg(bed.schedA, bed.indexB).withNotionals(bed.NB);
    Leg fixedB = FixedRateLeg(bed.schedB).withNotionals(bed.NB).withCouponRates(bed.R, bed.dc);
    auto eurOption = QuantLib::ext::make_shared<MultiLegOption>(
        std::vector<Leg>{fixedB, floatB}, std::vector<bool>{true, false},
        std::vector<Currency>{EURCurrency(), EURCurrency()}, QuantLib::ext::make_shared<BermudanExercise>(bed.exDates),
        Settlement::Physical, Settlement::PhysicalOTC, bed.exDates);
    eurOption->setPricingEngine(QuantLib::ext::make_shared<FmmLsmMultiLegOptionEngine>(bed.model, bed.grid, cfg));
    const Real single = eurOption->NPV();
    BOOST_TEST_MESSAGE("reduced engine " << reduced << " EUR vs single-currency EUR Bermudan " << single << " (diff "
                                         << reduced - single << ", LSM s.e. "
                                         << option->result<Real>("fmmLsmLowerBoundStdError") << ")");
    BOOST_CHECK_SMALL(reduced - single, 1e-6 * bed.NB);
    // NPV in the USD reporting currency: the EUR value times USD per EUR
    option->setPricingEngine(QuantLib::ext::make_shared<FmmLsmXccyReducedEngine>(bed.inputs, USDCurrency(), gridProvider,
                                                                                  modelProvider, cfg));
    BOOST_CHECK_SMALL(option->NPV() - reduced / bed.fxBperA->value(), 1e-8 * bed.NB);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
