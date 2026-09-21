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

#include <qle/models/fmmexposure.hpp>

#include <ql/currencies/europe.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <cmath>

using namespace QuantLib;
using namespace QuantExt;

BOOST_FIXTURE_TEST_SUITE(QuantExtTestSuite, qle::test::TopLevelFixture)

BOOST_AUTO_TEST_SUITE(FmmExposureTest)

namespace {

// quarterly grid to 5y on a flat 3% curve, one factor (the exposure engine is model-agnostic)
struct ExposureBed {
    ExposureBed(const Size factors = 1, const Real rhoInf = 1.0, const Real beta = 0.0) {
        Settings::instance().evaluationDate() = Date(19, September, 2026);
        curve = Handle<YieldTermStructure>(
            QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.03, Actual365Fixed()));
        const Size M = 20;
        Array rateTimes(M + 1);
        for (Size k = 0; k <= M; ++k)
            rateTimes[k] = 0.25 * static_cast<Real>(k);
        Array shifts(M, 1.0 / 0.25);
        std::vector<Array> volLevels(M, Array(1, 0.0025));
        parametrization = QuantLib::ext::make_shared<FmmParametrization>(
            EURCurrency(), curve, rateTimes, shifts, Array(), volLevels,
            FmmParametrization::LocalVolType::DisplacedDiffusion, rhoInf, beta, factors);
        model = QuantLib::ext::make_shared<ForwardMarketModel>(parametrization);
    }
    // payer swap flows (receive float quarterly, pay fixed K annually) over grid indices (from, to]
    FmmCallableInstrument swapFlows(const Size from, const Size to, const Real K) const {
        const auto& p = *parametrization;
        FmmCallableInstrument inst;
        inst.style = FmmCallableInstrument::Style::Cancel;
        inst.lastFlowIdx = to;
        inst.fixedFlows.assign(p.numberOfRates() + 1, 0.0);
        inst.floatWeights.assign(p.numberOfRates() + 1, 0.0);
        for (Size j = from + 1; j <= to; ++j)
            inst.floatWeights[j] = 1.0;
        for (Size c = from + 4; c <= to; c += 4)
            inst.fixedFlows[c] = -K * (p.rateTime(c) - p.rateTime(c - 4));
        return inst;
    }
    std::vector<Time> gridTimes(const Size from, const Size to) const {
        std::vector<Time> t;
        for (Size j = from; j <= to; ++j)
            t.push_back(parametrization->rateTime(j));
        return t;
    }
    Handle<YieldTermStructure> curve;
    QuantLib::ext::shared_ptr<FmmParametrization> parametrization;
    QuantLib::ext::shared_ptr<ForwardMarketModel> model;
};

FmmExposureConfig smallConfig() {
    FmmExposureConfig c;
    c.trainingPaths = 4096;
    c.trainingSeed = 4242;
    c.exposurePaths = 2000;
    c.exposureSeed = 424242;
    return c;
}

// a policy that exercises at the first in-the-money right (continuation estimate -1e12)
FmmLsmPolicy forcedPolicy(const Size nRights, const Size nb) {
    FmmLsmPolicy pol;
    pol.coefficients.assign(nRights, Array(nb, 0.0));
    for (auto& c : pol.coefficients)
        c[0] = -1e12;
    return pol;
}

} // namespace

BOOST_AUTO_TEST_CASE(testVanillaSwapExactValuesAndMartingale) {
    BOOST_TEST_MESSAGE("Testing the exposure engine on a vanilla swap: cube values equal the exact curve values on "
                       "every path and date (on and off the tenor grid), martingale identity, determinism...");
    ExposureBed bed;
    FmmExposureTrade swap;
    swap.id = "swap";
    swap.instrument = bed.swapFlows(4, 20, 0.03);
    std::vector<Time> times;
    for (Size k = 1; k <= 20; ++k) {
        times.push_back(0.25 * k - 0.125); // mid-period (off grid)
        times.push_back(0.25 * k);         // grid dates
    }
    const auto config = smallConfig();
    FmmExposureEngine engine(bed.model, {swap}, times, config);
    const auto res = engine.calculate();
    BOOST_REQUIRE_EQUAL(res.deflatedValues.size(), Size(1));
    BOOST_REQUIRE_EQUAL(res.times.size(), times.size());
    // exact values: regenerate the exposure paths (same user times = union of grid and exposure
    // times, same seed) and mark the flows on each state
    std::vector<Time> user;
    for (Size k = 1; k <= 20; ++k) {
        user.push_back(0.25 * k - 0.125);
        user.push_back(0.25 * k);
    }
    FmmPathGenerator gen(bed.model, user, config.exposureSequence, config.exposureSeed);
    Real maxErr = 0.0;
    for (Size n = 0; n < config.exposurePaths; ++n) {
        const auto path = gen.next();
        for (Size d = 0; d < times.size(); ++d) {
            const auto& st = path.states[d];
            const Real bank = bed.model->bankAccount(st);
            const Real exact = FmmExposureEngine::markValue(*bed.model, swap.instrument, st, 0, 20, false);
            maxErr = std::max(maxErr, std::fabs(res.deflatedValues[0][d * config.exposurePaths + n] * bank - exact));
            BOOST_CHECK_CLOSE(res.numeraire[d * config.exposurePaths + n], bank, 1e-10);
        }
    }
    BOOST_TEST_MESSAGE("max |cube value - exact| " << maxErr);
    BOOST_CHECK_SMALL(maxErr, 1e-12);
    // martingale identity per date: E[V(t)/B(t)] + past flows = V(0) (the curve value)
    const Real v0 = res.t0CurveValues[0];
    BOOST_CHECK_CLOSE(res.t0Values[0], v0, 100.0 * 3.0 * res.t0ValuesSe[0] / std::fabs(v0) + 1e-9);
    Size worstDate = 0;
    Real worstZ = 0.0;
    for (Size d = 0; d < times.size(); ++d) {
        const auto& diag = res.diagnostics[0][d];
        BOOST_CHECK_EQUAL(diag.aliveFraction, 0.0);
        const Real z = std::fabs(diag.martingaleCheck - v0) / diag.martingaleCheckSe;
        if (z > worstZ) {
            worstZ = z;
            worstDate = d;
        }
    }
    BOOST_TEST_MESSAGE("swap V(0) " << v0 << ", MC " << res.t0Values[0] << " +/- " << res.t0ValuesSe[0]
                                    << ", worst martingale z " << worstZ << " at t " << times[worstDate]);
    BOOST_CHECK(worstZ < 4.0);
    // values after the last flow are zero, before the first flow the whole swap
    const Size lastDate = times.size() - 1;
    for (Size n = 0; n < config.exposurePaths; ++n)
        BOOST_CHECK_SMALL(res.deflatedValues[0][lastDate * config.exposurePaths + n], 1e-14);
    // determinism
    const auto res2 = FmmExposureEngine(bed.model, {swap}, times, config).calculate();
    BOOST_CHECK_EQUAL(res2.deflatedValues[0][5 * config.exposurePaths + 7], res.deflatedValues[0][5 * config.exposurePaths + 7]);
    BOOST_CHECK_EQUAL(res2.numeraire[9 * config.exposurePaths + 3], res.numeraire[9 * config.exposurePaths + 3]);
}

BOOST_AUTO_TEST_CASE(testForcedExerciseEventRules) {
    BOOST_TEST_MESSAGE("Testing the contract-event rules: physically settled swaption exercised at the first "
                       "in-the-money right becomes the swap (with and without a notice period), a cancelled "
                       "swap dies after settlement, fees are marked at settlement...");
    ExposureBed bed;
    const auto config = smallConfig();
    FmmExposureTrade swap;
    swap.id = "swap";
    swap.instrument = bed.swapFlows(4, 20, 0.03);
    // Enter: exercise at 1y (right 0, notice = settle = 4), 2y (8), 3y (12) into the swap to 5y
    FmmExposureTrade berm;
    berm.id = "berm";
    berm.instrument = bed.swapFlows(4, 20, 0.03);
    berm.instrument.style = FmmCallableInstrument::Style::Enter;
    berm.instrument.rights = {{4, 4, 0.0}, {8, 8, 0.0}, {12, 12, 0.0}};
    // Enter with a one-period notice: decide at 4, enter at 5 (flows after index 5), fee 0.001 at 5
    FmmExposureTrade bermNotice = berm;
    bermNotice.id = "bermNotice";
    bermNotice.instrument.rights = {{4, 5, 0.001}, {8, 9, 0.001}, {12, 13, 0.001}};
    // Cancel: the swap holder may cancel at 4, 8, 12 paying a fee of 0.002 at settlement
    FmmExposureTrade cancel;
    cancel.id = "cancel";
    cancel.instrument = bed.swapFlows(4, 20, 0.03);
    cancel.instrument.rights = {{4, 4, -0.002}, {8, 8, -0.002}, {12, 12, -0.002}};
    FmmLsmPricer helper(bed.model, berm.instrument);
    const Size nb = helper.basis(Array(3, 0.0)).size();
    berm.policy = forcedPolicy(3, nb);
    bermNotice.policy = forcedPolicy(3, nb);
    cancel.policy = forcedPolicy(3, nb);
    const std::vector<Time> times = bed.gridTimes(1, 20);
    FmmExposureEngine engine(bed.model, {swap, berm, bermNotice, cancel}, times, config);
    const auto res = engine.calculate();
    const Size paths = config.exposurePaths;
    // regenerate the exposure paths (user times = the grid times, as in the engine) and classify
    FmmPathGenerator gen(bed.model, times, config.exposureSequence, config.exposureSeed);
    FmmLsmPricer bermHelper(bed.model, berm.instrument), noticeHelper(bed.model, bermNotice.instrument),
        cancelHelper(bed.model, cancel.instrument);
    Size exercisedAt0 = 0, checkedNotice = 0, checkedCancel = 0;
    Real maxErr = 0.0, maxErrNotice = 0.0, maxErrCancel = 0.0;
    for (Size n = 0; n < paths; ++n) {
        const auto path = gen.next();
        auto itm = [&](const FmmLsmPricer& h, const Size r) {
            const auto& rt = h.policy().coefficients.empty() ? berm.instrument.rights[r] : berm.instrument.rights[r];
            (void)rt;
            const Size notice = r == 0 ? 4 : r == 1 ? 8 : 12;
            const auto& st = path.states[notice - 1];
            return h.regressorsAt(st, r, bed.model->bankAccount(st))[2] > 0.0;
        };
        // Cancel: cancelled at T_4 (fee paid at T_4, excluded on that date): zero from T_4 on
        if (itm(cancelHelper, 0)) {
            ++checkedCancel;
            for (Size d = 3; d < times.size(); ++d)
                maxErrCancel = std::max(maxErrCancel, std::fabs(res.deflatedValues[3][d * paths + n]));
        }
        const bool ex0 = itm(bermHelper, 0);
        if (!ex0)
            continue;
        ++exercisedAt0;
        // Enter without notice: from the exercise date on, the swaption IS the swap
        for (Size d = 3; d < times.size(); ++d) { // d = 3 is t = 1y = T_4
            const Real vb = res.deflatedValues[1][d * paths + n], vs = res.deflatedValues[0][d * paths + n];
            maxErr = std::max(maxErr, std::fabs(vb - vs));
        }
        // Enter with notice: decided at T_4 (in the money there too? its intrinsic includes the
        // fee, check its own decision), swap delivered at T_5: at T_4 the committed value is the
        // flows after index 5 plus the fee at T_5, from T_5 on the swap flows after index 5
        if (itm(noticeHelper, 0)) {
            ++checkedNotice;
            const auto& st4 = path.states[3];
            const Real bank4 = bed.model->bankAccount(st4);
            const Real expected4 = FmmExposureEngine::markValue(*bed.model, bermNotice.instrument, st4, 5, 20, false) +
                                   FmmExposureEngine::feeValue(*bed.model, bermNotice.instrument, st4, 0, false);
            maxErrNotice = std::max(maxErrNotice, std::fabs(res.deflatedValues[2][3 * paths + n] * bank4 - expected4));
            for (Size d = 4; d < times.size(); ++d) { // T_5 on: the delivered swap, flows after index 5
                const auto& st = path.states[d];
                const Real bank = bed.model->bankAccount(st);
                const Real expected = FmmExposureEngine::markValue(*bed.model, bermNotice.instrument, st, 5, 20, false);
                maxErrNotice = std::max(maxErrNotice, std::fabs(res.deflatedValues[2][d * paths + n] * bank - expected));
            }
        }
    }
    BOOST_TEST_MESSAGE("paths exercised at the first right " << exercisedAt0 << " of " << paths << ": max |swaption - swap| "
                                                              << maxErr << "; notice variant checked on " << checkedNotice
                                                              << " paths, max error " << maxErrNotice
                                                              << "; cancel checked on " << checkedCancel
                                                              << " paths, max |value| " << maxErrCancel);
    BOOST_CHECK(exercisedAt0 > 100);
    BOOST_CHECK_SMALL(maxErr, 1e-12);
    BOOST_CHECK(checkedNotice > 100);
    BOOST_CHECK_SMALL(maxErrNotice, 1e-12);
    BOOST_CHECK(checkedCancel > 100);
    BOOST_CHECK_SMALL(maxErrCancel, 1e-12);
    // diagnostics: at T_4 nobody is alive for the no-notice swaption if exercised, else alive; the
    // exercised fraction is monotone in time and the notice variant is 'noticed' at T_4 only
    const auto& dg = res.diagnostics;
    BOOST_CHECK_CLOSE(dg[1][3].exercisedFraction, static_cast<Real>(exercisedAt0) / paths, 1e-9);
    BOOST_CHECK(dg[2][3].noticedFraction > 0.0);
    BOOST_CHECK_EQUAL(dg[2][4].noticedFraction, 0.0);
    for (Size d = 1; d < times.size(); ++d)
        BOOST_CHECK(dg[1][d].exercisedFraction >= dg[1][d - 1].exercisedFraction - 1e-12);
    // before the first notice the swaption is alive on every path (regression), after the last
    // notice nothing is pending
    BOOST_CHECK_EQUAL(dg[1][0].aliveFraction, 1.0);
    BOOST_CHECK_EQUAL(dg[1][12].aliveFraction, 0.0);
}

BOOST_AUTO_TEST_CASE(testCancelNeverIdentityAndTrainedPolicyMartingale) {
    BOOST_TEST_MESSAGE("Testing the conditional-value regression: a cancellable swap whose policy never cancels "
                       "must reproduce the vanilla swap (regression before the last notice, exact after), and the "
                       "martingale identity of a Bermudan swaption under a trained A4 policy...");
    ExposureBed bed;
    const auto config = smallConfig();
    FmmExposureTrade swap;
    swap.id = "swap";
    swap.instrument = bed.swapFlows(4, 20, 0.03);
    FmmExposureTrade never;
    never.id = "cancelNever";
    never.instrument = bed.swapFlows(4, 20, 0.03);
    never.instrument.rights = {{4, 4, 0.0}, {8, 8, 0.0}, {12, 12, 0.0}};
    never.policy.coefficients.assign(3, Array()); // empty = never exercise
    // Bermudan payer swaption with the A4-trained policy on the same model
    FmmExposureTrade berm;
    berm.id = "berm";
    berm.instrument = bed.swapFlows(4, 20, 0.03);
    berm.instrument.style = FmmCallableInstrument::Style::Enter;
    berm.instrument.rights = {{4, 4, 0.0}, {8, 8, 0.0}, {12, 12, 0.0}, {16, 16, 0.0}};
    FmmLsmConfig lsm;
    lsm.trainingPaths = 8192;
    lsm.valuationPaths = 8192;
    lsm.trainingSeed = 77;
    lsm.valuationSeed = 78;
    FmmLsmPricer pricer(bed.model, berm.instrument, lsm);
    const auto lsmRes = pricer.calculate();
    berm.policy = pricer.policy();
    std::vector<Time> times;
    for (Size k = 1; k <= 20; ++k) {
        times.push_back(0.25 * k - 0.125);
        times.push_back(0.25 * k);
    }
    FmmExposureEngine engine(bed.model, {swap, never, berm}, times, config);
    const auto res = engine.calculate();
    const Size paths = config.exposurePaths;
    // cancel-never vs swap: exact after the last notice (t >= 3y), regression before; the exact
    // underlying part is taken out of the regression, so what remains is the (here absent)
    // cancellation option fitted to path noise: the estimate deviates from the exact swap value by
    // the regression noise only, whose rms scales like 1 / sqrt(training paths) (checked with four
    // times the training paths); R2 of the residual regression is not meaningful and only reported
    auto regionStats = [&](const FmmExposureResult& r, Real& maxExact, Real& maxReg, Real& rmsReg, Real& worstR2) {
        maxExact = maxReg = rmsReg = 0.0;
        worstR2 = 1.0;
        Size count = 0;
        for (Size d = 0; d < times.size(); ++d) {
            for (Size n = 0; n < paths; ++n) {
                const Real diff = std::fabs(r.deflatedValues[1][d * paths + n] - r.deflatedValues[0][d * paths + n]);
                if (times[d] >= 3.0 - 1e-12)
                    maxExact = std::max(maxExact, diff);
                else {
                    maxReg = std::max(maxReg, diff);
                    rmsReg += diff * diff;
                    ++count;
                }
            }
            if (times[d] < 3.0 - 1e-12)
                worstR2 = std::min(worstR2, r.diagnostics[1][d].regressionR2);
        }
        rmsReg = std::sqrt(rmsReg / static_cast<Real>(std::max<Size>(count, 1)));
    };
    Real maxExact, maxReg, rmsReg, worstR2;
    regionStats(res, maxExact, maxReg, rmsReg, worstR2);
    auto config4 = config;
    config4.trainingPaths = 4 * config.trainingPaths;
    const auto res4 = FmmExposureEngine(bed.model, {swap, never, berm}, times, config4).calculate();
    Real maxExact4, maxReg4, rmsReg4, worstR24;
    regionStats(res4, maxExact4, maxReg4, rmsReg4, worstR24);
    BOOST_TEST_MESSAGE("cancel-never vs swap: exact region max |diff| " << maxExact << "; regression region rms "
                                                                        << rmsReg << ", max " << maxReg << " with "
                                                                        << config.trainingPaths << " training paths, rms "
                                                                        << rmsReg4 << ", max " << maxReg4 << " with "
                                                                        << config4.trainingPaths
                                                                        << " (deflated, unit notional); worst R2 "
                                                                        << worstR2);
    BOOST_CHECK_SMALL(maxExact, 1e-12);
    BOOST_CHECK_SMALL(maxExact4, 1e-12);
    BOOST_CHECK(rmsReg < 5e-3);
    BOOST_CHECK(maxReg < 5e-2);
    BOOST_CHECK(rmsReg4 < 0.75 * rmsReg); // noise scaling with the training paths
    // Bermudan: time-zero value on the exposure paths agrees with the LSM valuation, the martingale
    // identity holds at every date within the estimator noise, the out-of-sample bias is small
    BOOST_TEST_MESSAGE("Bermudan LSM lower bound " << lsmRes.lowerBound << " +/- " << lsmRes.lowerBoundSe
                                                  << ", exposure-path value " << res.t0Values[2] << " +/- "
                                                  << res.t0ValuesSe[2]);
    BOOST_CHECK(std::fabs(res.t0Values[2] - lsmRes.lowerBound) <
                3.0 * std::sqrt(res.t0ValuesSe[2] * res.t0ValuesSe[2] + lsmRes.lowerBoundSe * lsmRes.lowerBoundSe));
    Real worstZ = 0.0, worstBias = 0.0;
    Time worstT = 0.0;
    for (Size d = 0; d < times.size(); ++d) {
        const auto& dg = res.diagnostics[2][d];
        if (dg.martingaleCheckSe > 0.0) {
            const Real z = std::fabs(dg.martingaleCheck - res.t0Values[2]) / dg.martingaleCheckSe;
            if (z > worstZ) {
                worstZ = z;
                worstT = times[d];
            }
        }
        if (dg.outOfSampleBias != Null<Real>())
            worstBias = std::max(worstBias, std::fabs(dg.outOfSampleBias) - 3.0 * dg.outOfSampleBiasSe);
    }
    BOOST_TEST_MESSAGE("Bermudan martingale identity: worst z " << worstZ << " at t " << worstT
                                                                << ", worst (|bias| - 3 s.e.) " << worstBias);
    BOOST_CHECK(worstZ < 4.0);
    BOOST_CHECK(worstBias < 2e-4); // 2 bp of notional beyond three standard errors
    // exercise probabilities: after the last notice the alive fraction is zero and the exercised
    // fraction equals the LSM's total exercise probability up to noise
    Real lsmExercise = 0.0;
    for (const Real q : lsmRes.exerciseProbability)
        lsmExercise += q;
    const auto& last = res.diagnostics[2].back();
    BOOST_TEST_MESSAGE("exercised fraction at the end " << last.exercisedFraction << ", LSM total exercise probability "
                                                        << lsmExercise);
    BOOST_CHECK_EQUAL(last.aliveFraction, 0.0);
    BOOST_CHECK_SMALL(last.exercisedFraction - lsmExercise, 0.03);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
