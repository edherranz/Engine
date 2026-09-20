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

#include <ql/currencies/europe.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/integrals/gaussianquadratures.hpp>
#include <ql/math/solvers1d/brent.hpp>
#include <ql/math/statistics/incrementalstatistics.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <qle/instruments/callablebond.hpp>
#include <qle/models/fmmanalytics.hpp>
#include <qle/models/fmmcalibration.hpp>
#include <qle/models/fmmlsmpricer.hpp>
#include <qle/models/fmmparametrization.hpp>
#include <qle/models/forwardmarketmodel.hpp>
#include <qle/models/irlgm1fconstantparametrization.hpp>
#include <qle/models/irlgm1fpiecewiseconstantparametrization.hpp>
#include <qle/models/lgm.hpp>
#include <qle/pricingengines/analyticlgmswaptionengine.hpp>
#include <qle/pricingengines/numericlgmcallablebondengine.hpp>
#include <qle/pricingengines/numericlgmmultilegoptionengine.hpp>

#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/instruments/makevanillaswap.hpp>
#include <ql/instruments/swaption.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/models/shortrate/calibrationhelpers/swaptionhelper.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/time/schedule.hpp>

#include <cstdlib>

using namespace QuantLib;
using namespace QuantExt;

BOOST_FIXTURE_TEST_SUITE(QuantExtTestSuite, qle::test::TopLevelFixture)

BOOST_AUTO_TEST_SUITE(FmmModelTest)

namespace {

struct FmmTestBed {
    // quarterly grid to 5y on a flat 3% (continuous) curve
    FmmTestBed(const FmmParametrization::LocalVolType volType, const Real level, const Real shiftValue,
               const Size factors, const Real rhoInf = 0.6, const Real beta = 0.08) {
        Settings::instance().evaluationDate() = Date(19, September, 2026);
        curve = Handle<YieldTermStructure>(
            QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.03, Actual365Fixed()));
        const Size M = 20;
        Array rateTimes(M + 1);
        for (Size k = 0; k <= M; ++k)
            rateTimes[k] = 0.25 * static_cast<Real>(k);
        Array shifts(M);
        for (Size j = 0; j < M; ++j)
            shifts[j] = shiftValue < 0.0 ? 1.0 / 0.25 : shiftValue; // negative marker = 1/tau
        Array volTimes; // flat levels
        std::vector<Array> volLevels(M, Array(1, level));
        parametrization = QuantLib::ext::make_shared<FmmParametrization>(
            EURCurrency(), curve, rateTimes, shifts, volTimes, volLevels, volType, rhoInf, beta, factors);
        model = QuantLib::ext::make_shared<ForwardMarketModel>(parametrization);
    }

    std::vector<Time> allGridTimes() const {
        std::vector<Time> t;
        for (Size j = 1; j <= parametrization->numberOfRates(); ++j)
            t.push_back(parametrization->rateTime(j));
        return t;
    }

    Handle<YieldTermStructure> curve;
    QuantLib::ext::shared_ptr<FmmParametrization> parametrization;
    QuantLib::ext::shared_ptr<ForwardMarketModel> model;
};

// analytic caplet on period j using the parametrization's exact integrated variance
Real analyticCaplet(const FmmParametrization& p, const Handle<YieldTermStructure>& curve, const Size j, const Real K,
                    const bool backwardLooking) {
    const Real tau = p.tau(j);
    const Real F0 = (curve->discount(p.rateTime(j - 1)) / curve->discount(p.rateTime(j)) - 1.0) / tau;
    const Time Te = backwardLooking ? p.rateTime(j) : p.rateTime(j - 1);
    const Real var = p.integratedCovariance(j, j, 0.0, Te) / p.correlation(j, j);
    const Real df = curve->discount(p.rateTime(j));
    if (p.volType() == FmmParametrization::LocalVolType::Normal)
        return df * tau * bachelierBlackFormula(Option::Call, K, F0, std::sqrt(var));
    return df * tau * blackFormula(Option::Call, K, F0, std::sqrt(var), 1.0, p.shift(j));
}

} // namespace

BOOST_AUTO_TEST_CASE(testMartingaleGridAndOffGrid) {
    BOOST_TEST_MESSAGE("Testing FMM deflated bond martingale property on and off the tenor grid...");

    FmmTestBed bed(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, 3);
    const Size paths = 16384;

    // user times: all grid dates plus two off-grid times
    std::vector<Time> userTimes = bed.allGridTimes();
    const Time tOff1 = 1.1, tOff2 = 3.6;
    userTimes.push_back(tOff1);
    userTimes.push_back(tOff2);
    std::sort(userTimes.begin(), userTimes.end());

    FmmPathGenerator gen(bed.model, userTimes, MersenneTwister, 42);

    // deflated bonds at a lattice of (observation, maturity) pairs, incl. off-grid maturities
    struct Probe {
        Time obs, mat;
        IncrementalStatistics stats;
    };
    std::vector<Probe> probes;
    for (const Time obs : {0.5, 1.0, tOff1, 2.0, tOff2, 4.0}) {
        for (const Time mat : {1.25, 1.6, 3.0, 4.3, 5.0}) {
            if (mat > obs + 0.01) {
                Probe p;
                p.obs = obs;
                p.mat = mat;
                probes.push_back(p);
            }
        }
    }

    for (Size n = 0; n < paths; ++n) {
        const auto path = gen.next();
        for (auto& probe : probes) {
            const Size idx = std::lower_bound(userTimes.begin(), userTimes.end(), probe.obs - 1e-10) -
                             userTimes.begin();
            const auto& state = path.states[idx];
            probe.stats.add(bed.model->discountBond(state, probe.mat) / bed.model->bankAccount(state));
        }
    }

    for (auto& probe : probes) {
        const Real expected = bed.curve->discount(probe.mat);
        const Real err = std::max(probe.stats.errorEstimate(), 1e-12);
        const Real diff = probe.stats.mean() - expected;
        BOOST_TEST_MESSAGE("E[P(" << probe.obs << "," << probe.mat << ")/B] = " << probe.stats.mean() << " vs "
                                  << expected << " diff " << diff << " s.e. " << err);
        BOOST_CHECK_MESSAGE(std::fabs(diff) < 3.0 * err + 5e-5,
                            "martingale violation at obs " << probe.obs << " mat " << probe.mat << ": diff " << diff
                                                           << " vs 3 s.e. " << 3.0 * err);
    }
}

BOOST_AUTO_TEST_CASE(testCapletsAgainstClosedForm) {
    BOOST_TEST_MESSAGE("Testing FMM Monte Carlo caplets against closed forms (DD and normal)...");

    struct Config {
        FmmParametrization::LocalVolType type;
        Real level, shift;
        const char* label;
    };
    const std::vector<Config> configs = {
        {FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, "DD shift 1/tau"},
        {FmmParametrization::LocalVolType::DisplacedDiffusion, 0.15, 0.02, "DD shift 2%"},
        {FmmParametrization::LocalVolType::Normal, 0.008, 0.0, "normal"}};

    const Size paths = 16384;
    for (const auto& cfg : configs) {
        FmmTestBed bed(cfg.type, cfg.level, cfg.shift, 3);
        std::vector<Time> userTimes = bed.allGridTimes();
        FmmPathGenerator gen(bed.model, userTimes, MersenneTwister, 7);

        const std::vector<Size> testRates = {4, 12, 20};
        const Real F0 = 0.03; // approximate atm level on the flat curve
        const std::vector<Real> strikes = {F0 - 0.01, F0, F0 + 0.01};

        std::map<std::pair<Size, Size>, IncrementalStatistics> statsB, statsF; // (rate, strikeIdx)
        for (Size n = 0; n < paths; ++n) {
            const auto path = gen.next();
            for (const Size j : testRates) {
                const auto& stateEnd = path.states[j - 1];        // at T_j
                const auto& stateFix = j >= 2 ? path.states[j - 2] : path.states[0]; // at T_{j-1} (j>=2)
                const Real bEnd = bed.model->bankAccount(stateEnd);
                for (Size ki = 0; ki < strikes.size(); ++ki) {
                    const Real tau = bed.parametrization->tau(j);
                    statsB[{j, ki}].add(tau * std::max(stateEnd.R[j - 1] - strikes[ki], 0.0) / bEnd);
                    if (j >= 2)
                        statsF[{j, ki}].add(tau * std::max(stateFix.R[j - 1] - strikes[ki], 0.0) / bEnd);
                }
            }
        }

        for (const Size j : testRates) {
            for (Size ki = 0; ki < 3; ++ki) {
                const Real anB = analyticCaplet(*bed.parametrization, bed.curve, j, strikes[ki], true);
                const Real mcB = statsB[{j, ki}].mean();
                const Real seB = statsB[{j, ki}].errorEstimate();
                BOOST_TEST_MESSAGE(cfg.label << " backward caplet j=" << j << " K=" << strikes[ki] << ": mc " << mcB
                                             << " analytic " << anB << " diff " << mcB - anB << " s.e. " << seB);
                BOOST_CHECK_MESSAGE(std::fabs(mcB - anB) < 3.0 * seB + 2e-6,
                                    cfg.label << " backward caplet j=" << j << " K=" << strikes[ki] << " diff "
                                              << mcB - anB << " exceeds 3 s.e. " << 3.0 * seB);
                if (j >= 2) {
                    const Real anF = analyticCaplet(*bed.parametrization, bed.curve, j, strikes[ki], false);
                    const Real mcF = statsF[{j, ki}].mean();
                    const Real seF = statsF[{j, ki}].errorEstimate();
                    BOOST_CHECK_MESSAGE(std::fabs(mcF - anF) < 3.0 * seF + 2e-6,
                                        cfg.label << " forward caplet j=" << j << " K=" << strikes[ki] << " diff "
                                                  << mcF - anF << " exceeds 3 s.e. " << 3.0 * seF);
                    // model-free ordering: backward >= forward
                    BOOST_CHECK_MESSAGE(anB > anF - 1e-12, "backward caplet cheaper than forward");
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(testMeasureInvariance) {
    BOOST_TEST_MESSAGE("Testing FMM European price invariance between spot and terminal measures...");

    const Size j = 8; // caplet on [T_7, T_8]
    const Real K = 0.03;
    const Size paths = 16384;

    FmmTestBed bedSpot(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, 3);
    std::vector<Time> userTimes = bedSpot.allGridTimes();

    // spot measure
    FmmPathGenerator genSpot(bedSpot.model, userTimes, MersenneTwister, 99);
    IncrementalStatistics statsSpot;
    for (Size n = 0; n < paths; ++n) {
        const auto path = genSpot.next();
        const auto& s = path.states[j - 1];
        statsSpot.add(bedSpot.parametrization->tau(j) * std::max(s.R[j - 1] - K, 0.0) /
                      bedSpot.model->bankAccount(s));
    }

    // terminal measure Q^{T_M}
    FmmTestBed bedTerm(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, 3);
    const Size M = bedTerm.parametrization->numberOfRates();
    auto modelTerm = QuantLib::ext::make_shared<ForwardMarketModel>(
        bedTerm.parametrization, ForwardMarketModel::Measure::TerminalForward, M);
    FmmPathGenerator genTerm(modelTerm, userTimes, MersenneTwister, 99);
    IncrementalStatistics statsTerm;
    const Real PTM = bedTerm.curve->discount(bedTerm.parametrization->rateTime(M));
    for (Size n = 0; n < paths; ++n) {
        const auto path = genTerm.next();
        const auto& s = path.states[j - 1];
        statsTerm.add(PTM * bedTerm.parametrization->tau(j) * std::max(s.R[j - 1] - K, 0.0) /
                      modelTerm->numeraire(s));
    }

    const Real diff = statsSpot.mean() - statsTerm.mean();
    const Real err = std::sqrt(statsSpot.errorEstimate() * statsSpot.errorEstimate() +
                               statsTerm.errorEstimate() * statsTerm.errorEstimate());
    BOOST_TEST_MESSAGE("spot " << statsSpot.mean() << " terminal " << statsTerm.mean() << " diff " << diff << " s.e. "
                               << err);
    BOOST_CHECK_MESSAGE(std::fabs(diff) < 3.0 * err + 2e-6,
                        "measure invariance violated: diff " << diff << " vs 3 s.e. " << 3.0 * err);
}

BOOST_AUTO_TEST_CASE(testLgmReplication) {
    BOOST_TEST_MESSAGE("Testing FMM replication mode against ORE's 1F LGM analytics...");

    Settings::instance().evaluationDate() = Date(19, September, 2026);
    Handle<YieldTermStructure> curve(
        QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.03, Actual365Fixed()));
    const Real alpha = 0.01, kappa = 0.01;
    auto lgm = QuantLib::ext::make_shared<IrLgm1fConstantParametrization>(EURCurrency(), curve, alpha, kappa);

    const Size M = 20;
    Array rateTimes(M + 1);
    for (Size k = 0; k <= M; ++k)
        rateTimes[k] = 0.25 * static_cast<Real>(k);
    auto fmmParam = QuantLib::ext::make_shared<FmmParametrization>(EURCurrency(), curve, rateTimes, lgm);
    auto model = QuantLib::ext::make_shared<ForwardMarketModel>(fmmParam);

    const Time Te = 2.0, Tm = 5.0;
    std::vector<Time> userTimes = {Te};
    const Size paths = 16384;
    FmmPathGenerator gen(model, userTimes, MersenneTwister, 4711);

    // LGM references from ORE's own analytics: the built-in zero-bond option formula
    // (LGM::discountBondOption with t = S = Te) and the ln-bond variance zeta(Te) (H(Tm)-H(Te))^2
    LinearGaussMarkovModel lgmModel(lgm);
    const Real zetaTe = lgm->zeta(Te);
    const Real K = curve->discount(Tm) / curve->discount(Te); // atm forward bond price
    const Real lgmPut = lgmModel.discountBondOption(Option::Put, K, Te, Te, Tm);

    IncrementalStatistics mcPut, mcLogBond;
    for (Size n = 0; n < paths; ++n) {
        const auto path = gen.next();
        const auto& s = path.states[0];
        const Real p = model->discountBond(s, Tm);
        mcPut.add(std::max(K - p, 0.0) / model->bankAccount(s));
        mcLogBond.add(std::log(p));
    }

    const Real diffPut = mcPut.mean() - lgmPut;
    BOOST_TEST_MESSAGE("zero-bond put: fmm-mc " << mcPut.mean() << " lgm " << lgmPut << " diff " << diffPut << " s.e. "
                                                << mcPut.errorEstimate());
    BOOST_CHECK_MESSAGE(std::fabs(diffPut) < 3.0 * mcPut.errorEstimate() + 2e-6,
                        "replication zero-bond put mismatch: " << diffPut);

    // moment check for ln P(Te, Tm): variance must equal (H(Tm)-H(Te))^2 zeta(Te)
    const Real lgmVar = std::pow(lgm->H(Tm) - lgm->H(Te), 2.0) * zetaTe;
    const Real mcVar = mcLogBond.variance();
    BOOST_TEST_MESSAGE("ln P variance: fmm-mc " << mcVar << " lgm " << lgmVar << " rel diff "
                                                << (mcVar - lgmVar) / lgmVar);
    BOOST_CHECK_MESSAGE(std::fabs(mcVar - lgmVar) / lgmVar < 0.03,
                        "replication ln-bond variance mismatch: " << mcVar << " vs " << lgmVar);
}

BOOST_AUTO_TEST_CASE(testDeterminism) {
    BOOST_TEST_MESSAGE("Testing FMM path generation determinism (same seed, bitwise equality)...");

    FmmTestBed bed(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, 3);
    std::vector<Time> userTimes = bed.allGridTimes();
    FmmPathGenerator gen1(bed.model, userTimes, SobolBrownianBridge, 12345);
    FmmPathGenerator gen2(bed.model, userTimes, SobolBrownianBridge, 12345);
    for (Size n = 0; n < 64; ++n) {
        const auto p1 = gen1.next();
        const auto p2 = gen2.next();
        for (Size i = 0; i < p1.states.size(); ++i) {
            for (Size j = 0; j < p1.states[i].R.size(); ++j)
                BOOST_REQUIRE(p1.states[i].R[j] == p2.states[i].R[j]);
            BOOST_REQUIRE(p1.states[i].x == p2.states[i].x);
        }
    }
    BOOST_TEST_MESSAGE("64 paths bitwise identical");
}

BOOST_AUTO_TEST_CASE(testStepBias) {
    BOOST_TEST_MESSAGE("Testing FMM discretisation bias: quarterly steps vs monthly sub-steps...");

    // strongly skewed configuration to expose drift discretisation error
    FmmTestBed bed(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.15, 0.02, 3);
    const Size j = 20;
    const Real K = 0.03;
    const Real analytic = analyticCaplet(*bed.parametrization, bed.curve, j, K, true);
    const Size paths = 16384;

    for (const bool fine : {false, true}) {
        std::vector<Time> userTimes;
        if (fine) {
            for (Size k = 1; k <= 60; ++k)
                userTimes.push_back(5.0 * static_cast<Real>(k) / 60.0);
        } else {
            userTimes = bed.allGridTimes();
        }
        FmmPathGenerator gen(bed.model, userTimes, MersenneTwister, 31);
        IncrementalStatistics stats;
        for (Size n = 0; n < paths; ++n) {
            const auto path = gen.next();
            const auto& s = path.states.back();
            stats.add(bed.parametrization->tau(j) * std::max(s.R[j - 1] - K, 0.0) /
                      bed.model->bankAccount(s));
        }
        const Real diff = stats.mean() - analytic;
        BOOST_TEST_MESSAGE((fine ? "monthly " : "quarterly ") << "steps: mc " << stats.mean() << " analytic "
                                                              << analytic << " diff " << diff << " s.e. "
                                                              << stats.errorEstimate());
        BOOST_CHECK_MESSAGE(std::fabs(diff) < 3.0 * stats.errorEstimate() + 5e-6,
                            "discretisation bias with " << (fine ? "monthly" : "quarterly") << " steps: " << diff);
    }
}

BOOST_AUTO_TEST_CASE(testCompletionContract) {
    BOOST_TEST_MESSAGE("Testing FMM completion contract: curve recovery, P(t,t)=1, grid identity, "
                       "continuity, positivity (FMM_SPEC.md section 4.3)...");

    // zero volatility: every path is deterministic and the completion must reproduce the initial
    // curve exactly through all stubs, snapshots and telescoping
    FmmTestBed bed(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0, 0.02, 2);

    // initial state: exact curve recovery
    const auto s0 = bed.model->initialState();
    for (const Time T : {0.1, 0.25, 0.4, 1.0, 2.3, 3.75, 5.0}) {
        const Real p = bed.model->discountBond(s0, T);
        const Real p0 = bed.curve->discount(T);
        BOOST_CHECK_MESSAGE(std::fabs(p - p0) < 1e-12,
                            "initial curve recovery failed at T=" << T << ": " << p << " vs " << p0);
    }
    BOOST_CHECK_SMALL(bed.model->bankAccount(s0) - 1.0, 1e-14);

    // deterministic path: states straddling tenor dates and off-grid times
    std::vector<Time> userTimes = {0.2, 0.25, 0.3, 1.0, 1.1, 2.499, 2.5, 2.501, 4.0, 4.9};
    FmmPathGenerator gen(bed.model, userTimes, MersenneTwister, 1);
    const auto path = gen.next();
    for (Size i = 0; i < userTimes.size(); ++i) {
        const auto& st = path.states[i];
        const Time t = userTimes[i];
        // P(t,t) = 1
        BOOST_CHECK_MESSAGE(std::fabs(bed.model->discountBond(st, t) - 1.0) < 1e-12,
                            "P(t,t) != 1 at t=" << t);
        // curve recovery along the zero-vol path: P(t,T) = P0(T)/P0(t), B(t) = 1/P0(t)
        for (const Time T : {t + 0.05, t + 0.8, 5.0}) {
            if (T > 5.0 + 1e-12)
                continue;
            const Real p = bed.model->discountBond(st, T);
            const Real expected = bed.curve->discount(T) / bed.curve->discount(t);
            BOOST_CHECK_MESSAGE(std::fabs(p - expected) < 1e-11,
                                "zero-vol curve recovery failed at t=" << t << " T=" << T << ": " << p << " vs "
                                                                       << expected);
            BOOST_CHECK_MESSAGE(p > 0.0, "non-positive bond at t=" << t << " T=" << T);
        }
        const Real b = bed.model->bankAccount(st);
        BOOST_CHECK_MESSAGE(std::fabs(b - 1.0 / bed.curve->discount(t)) < 1e-11,
                            "zero-vol bank account mismatch at t=" << t);
    }

    // continuity across the tenor date 2.5: compare P(t,5) just before and just after
    const Size iBefore = 5, iAt = 6, iAfter = 7; // 2.499, 2.5, 2.501
    const Real pBefore = bed.model->discountBond(path.states[iBefore], 5.0);
    const Real pAt = bed.model->discountBond(path.states[iAt], 5.0);
    const Real pAfter = bed.model->discountBond(path.states[iAfter], 5.0);
    BOOST_CHECK_MESSAGE(std::fabs(pBefore - pAt) < 5e-4 && std::fabs(pAt - pAfter) < 5e-4,
                        "discontinuity at tenor date: " << pBefore << " / " << pAt << " / " << pAfter);

    // grid identity at a grid-time state with NON-zero vol: telescoping product equals completion
    FmmTestBed bedVol(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, 2);
    std::vector<Time> gridTimes = bedVol.allGridTimes();
    FmmPathGenerator genVol(bedVol.model, gridTimes, MersenneTwister, 2);
    const auto pathVol = genVol.next();
    const auto& s8 = pathVol.states[7]; // t = T_8 = 2.0
    for (Size m = 9; m <= 20; ++m) {
        Real prod = 1.0;
        for (Size i = 9; i <= m; ++i)
            prod /= 1.0 + bedVol.parametrization->tau(i) * s8.R[i - 1];
        const Real p = bedVol.model->discountBond(s8, bedVol.parametrization->rateTime(m));
        BOOST_CHECK_MESSAGE(std::fabs(p - prod) < 1e-11,
                            "grid identity violated at T_m, m=" << m << ": " << p << " vs " << prod);
        BOOST_CHECK_MESSAGE(p > 0.0, "non-positive grid bond at m=" << m);
    }
}

BOOST_AUTO_TEST_CASE(testBiasTableCRN) {
    BOOST_TEST_MESSAGE("Paired common-random-numbers bias table (FMM_SPEC.md section 9.7): "
                       "sub-steps per quarter in {1,2,3,6}, backward caplet j=20...");

    // deliberately skewed configuration (state-dependent drift factors) to expose stepping bias
    FmmTestBed bed(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.15, 0.02, 3);
    const Size M = bed.parametrization->numberOfRates();
    const Size j = 20;
    const Real K = 0.03;
    const Real tau = bed.parametrization->tau(j);
    const Real analytic = analyticCaplet(*bed.parametrization, bed.curve, j, K, true);

    // fine grid: 6 sub-steps per quarter (h = 1/24); coarser configs aggregate {6,3,2,1} fine steps
    const Size fineSub = 6, nFine = 20 * fineSub;
    std::vector<Time> fineTimes(nFine);
    for (Size i = 0; i < nFine; ++i)
        fineTimes[i] = 5.0 * static_cast<Real>(i + 1) / static_cast<Real>(nFine);
    std::vector<ForwardMarketModel::StepData> fineSteps;
    Time prev = 0.0;
    for (const Time t : fineTimes) {
        fineSteps.push_back(bed.model->makeStep(prev, t));
        prev = t;
    }
    const std::vector<Size> aggregate = {6, 3, 2, 1}; // fine steps per coarse step
    std::vector<std::vector<ForwardMarketModel::StepData>> cfgSteps;
    for (const Size a : aggregate) {
        std::vector<ForwardMarketModel::StepData> steps;
        prev = 0.0;
        for (Size i = a; i <= nFine; i += a) {
            steps.push_back(bed.model->makeStep(prev, fineTimes[i - 1]));
            prev = fineTimes[i - 1];
        }
        cfgSteps.push_back(steps);
    }

    const Size paths = 8192;
    auto variates = makeMultiPathVariateGenerator(MersenneTwister, M + 1, nFine, 20260919,
                                                  SobolBrownianGenerator::Steps, SobolRsg::JoeKuoD7);
    std::vector<IncrementalStatistics> price(aggregate.size());
    std::vector<IncrementalStatistics> pairedDiff(aggregate.size()); // vs finest config
    for (Size n = 0; n < paths; ++n) {
        const auto sample = variates->next();
        // correlated fine shocks
        std::vector<Array> vFine(nFine, Array(M + 1, 0.0));
        for (Size i = 0; i < nFine; ++i)
            for (Size a = 0; a <= M; ++a) {
                Real acc = 0.0;
                for (Size b = 0; b <= M; ++b)
                    acc += fineSteps[i].shockSqrt[a][b] * sample.value[i][b];
                vFine[i][a] = acc;
            }
        std::vector<Real> payoff(aggregate.size());
        for (Size c = 0; c < aggregate.size(); ++c) {
            const Size a = aggregate[c];
            ForwardMarketModel::State state = bed.model->initialState();
            Size fineIdx = 0;
            for (const auto& step : cfgSteps[c]) {
                Array v(M + 1, 0.0);
                for (Size k = 0; k < a; ++k, ++fineIdx)
                    for (Size q = 0; q <= M; ++q)
                        v[q] += vFine[fineIdx][q];
                bed.model->evolveWithCorrelatedShocks(state, step, v);
            }
            payoff[c] = tau * std::max(state.R[j - 1] - K, 0.0) / bed.model->bankAccount(state);
            price[c].add(payoff[c]);
        }
        for (Size c = 0; c + 1 < aggregate.size(); ++c)
            pairedDiff[c].add(payoff[c] - payoff.back());
    }

    BOOST_TEST_MESSAGE("BIAS-TABLE (analytic " << analytic << "):");
    BOOST_TEST_MESSAGE("substeps/quarter |    price    | price-analytic |   s.e.   | paired diff vs finest | s.e.(diff)");
    const std::vector<Size> subPerQuarter = {1, 2, 3, 6};
    for (Size c = 0; c < aggregate.size(); ++c) {
        const bool finest = c + 1 == aggregate.size();
        BOOST_TEST_MESSAGE(subPerQuarter[c]
                           << " | " << price[c].mean() << " | " << price[c].mean() - analytic << " | "
                           << price[c].errorEstimate() << " | " << (finest ? 0.0 : pairedDiff[c].mean()) << " | "
                           << (finest ? 0.0 : pairedDiff[c].errorEstimate()));
    }
    // finest config must agree with the closed form
    BOOST_CHECK_MESSAGE(std::fabs(price.back().mean() - analytic) < 3.0 * price.back().errorEstimate() + 2e-6,
                        "finest config bias " << price.back().mean() - analytic);
    // two or more sub-steps per quarter: paired bias below resolution (OQ4 evidence)
    for (Size c = 1; c + 1 < aggregate.size(); ++c)
        BOOST_CHECK_MESSAGE(std::fabs(pairedDiff[c].mean()) < 3.0 * pairedDiff[c].errorEstimate() + 2e-6,
                            subPerQuarter[c] << " sub-steps paired bias " << pairedDiff[c].mean() << " exceeds "
                                             << 3.0 * pairedDiff[c].errorEstimate());
}

BOOST_AUTO_TEST_CASE(testReplicationMomentsWithUncertainty) {
    BOOST_TEST_MESSAGE("Replication ln-bond variance with replication-based uncertainty, plus a "
                       "KS normality check (supplementary evidence)...");

    Settings::instance().evaluationDate() = Date(19, September, 2026);
    Handle<YieldTermStructure> curve(
        QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.03, Actual365Fixed()));
    auto lgm = QuantLib::ext::make_shared<IrLgm1fConstantParametrization>(EURCurrency(), curve, 0.01, 0.01);
    const Size M = 20;
    Array rateTimes(M + 1);
    for (Size k = 0; k <= M; ++k)
        rateTimes[k] = 0.25 * static_cast<Real>(k);
    auto param = QuantLib::ext::make_shared<FmmParametrization>(EURCurrency(), curve, rateTimes, lgm);
    auto model = QuantLib::ext::make_shared<ForwardMarketModel>(param);

    const Time Te = 2.0, Tm = 5.0;
    const Real lgmVar = std::pow(lgm->H(Tm) - lgm->H(Te), 2.0) * lgm->zeta(Te);

    const Size reps = 16, pathsPerRep = 4096;
    IncrementalStatistics varAcross;
    std::vector<Real> firstRepLogs;
    for (Size r = 0; r < reps; ++r) {
        FmmPathGenerator gen(model, {Te}, MersenneTwister, 1000 + r);
        IncrementalStatistics logs;
        for (Size n = 0; n < pathsPerRep; ++n) {
            const auto path = gen.next();
            const Real lp = std::log(model->discountBond(path.states[0], Tm));
            logs.add(lp);
            if (r == 0)
                firstRepLogs.push_back(lp);
        }
        varAcross.add(logs.variance());
    }
    const Real varMean = varAcross.mean();
    const Real varSe = varAcross.errorEstimate(); // s.e. of the mean across independent replications
    BOOST_TEST_MESSAGE("ln P(2,5) variance: fmm " << varMean << " +/- " << varSe << " (16 reps x 4096), lgm analytic "
                                                  << lgmVar << ", diff " << varMean - lgmVar);
    BOOST_CHECK_MESSAGE(std::fabs(varMean - lgmVar) < 3.0 * varSe,
                        "replication variance off by " << varMean - lgmVar << " vs 3 s.e. " << 3.0 * varSe);

    // supplementary: one-sample KS against the normal with the LGM variance (mean taken empirical;
    // the drift of ln P depends on the measure while variance and normality do not)
    std::sort(firstRepLogs.begin(), firstRepLogs.end());
    Real meanHat = 0.0;
    for (const Real v : firstRepLogs)
        meanHat += v;
    meanHat /= static_cast<Real>(firstRepLogs.size());
    CumulativeNormalDistribution N;
    Real d = 0.0;
    const Real n = static_cast<Real>(firstRepLogs.size());
    for (Size i = 0; i < firstRepLogs.size(); ++i) {
        const Real z = (firstRepLogs[i] - meanHat) / std::sqrt(lgmVar);
        const Real F = N(z);
        d = std::max(d, std::max(std::fabs(F - static_cast<Real>(i) / n),
                                 std::fabs(F - static_cast<Real>(i + 1) / n)));
    }
    BOOST_TEST_MESSAGE("KS statistic D_n = " << d << " (n = " << firstRepLogs.size()
                                             << "), 1% asymptotic threshold ~ " << 1.63 / std::sqrt(n));
    BOOST_CHECK_MESSAGE(d < 2.5 / std::sqrt(n), "KS distance " << d << " unexpectedly large");
}

BOOST_AUTO_TEST_CASE(testStochasticCompletionAndBoundary) {
    BOOST_TEST_MESSAGE("Stochastic completion identities and grid-boundary behaviour incl. T_M...");

    FmmTestBed bed(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, 2);
    const std::vector<Time> userTimes = {2.4999, 2.5, 2.5001, 3.7, 4.9999, 5.0};
    FmmPathGenerator gen(bed.model, userTimes, MersenneTwister, 5);
    const Size paths = 64;
    for (Size np = 0; np < paths; ++np) {
        const auto path = gen.next();
        for (Size i = 0; i < userTimes.size(); ++i) {
            const auto& st = path.states[i];
            const Time t = userTimes[i];
            BOOST_REQUIRE_MESSAGE(std::fabs(bed.model->discountBond(st, t) - 1.0) < 1e-10,
                                  "stochastic P(t,t) != 1 at t=" << t);
            for (const Time T : {4.0, 5.0}) {
                if (T < t)
                    continue;
                BOOST_REQUIRE_MESSAGE(bed.model->discountBond(st, T) > 0.0,
                                      "non-positive stochastic bond at t=" << t << " T=" << T);
            }
        }
        // continuity of ln P(., 4.0) across the tenor date 2.5 (eps = 1e-4 on both sides)
        const Real lpBefore = std::log(bed.model->discountBond(path.states[0], 4.0));
        const Real lpAfter = std::log(bed.model->discountBond(path.states[2], 4.0));
        BOOST_REQUIRE_MESSAGE(std::fabs(lpAfter - lpBefore) < 0.02,
                              "ln P jump across tenor date: " << lpAfter - lpBefore);
        // T_M boundary: state rests exactly at the last tenor date
        const auto& sM = path.states.back();
        BOOST_REQUIRE(std::fabs(sM.t - 5.0) < 1e-12);
        BOOST_REQUIRE_MESSAGE(std::fabs(bed.model->discountBond(sM, 5.0) - 1.0) < 1e-10, "P(T_M,T_M) != 1");
        Real prod = 1.0;
        for (Size jj = 1; jj <= 20; ++jj)
            prod *= 1.0 + bed.parametrization->tau(jj) * sM.R[jj - 1];
        const Real b = bed.model->bankAccount(sM);
        BOOST_REQUIRE_MESSAGE(std::fabs(b / prod - 1.0) < 1e-10,
                              "B(T_M) grid identity violated: " << b << " vs " << prod);
        BOOST_REQUIRE(bed.model->numeraire(sM) > 0.0);
    }
    BOOST_TEST_MESSAGE("64 stochastic paths: P(t,t)=1, positivity, tenor-date continuity, exact "
                       "B(T_M) grid identity all hold");
}

namespace {

// ---- test-local skew calibration machinery (library-grade calibrators are an A3 deliverable;
// ---- per house test convention these helpers are not called from the library) ----

// backward-caplet implied normal vol (annualised over [0, T_j]) for bucket j under DD (lambda, delta)
Real ddCapletNormalVol(const Real F0, const Real K, const Real lambda, const Real delta, const Real I,
                       const Time tte) {
    const Real stdDev = lambda * std::sqrt(I); // I = int_0^{T_j} g_j^2 du, exact
    const Real price = blackFormula(Option::Call, K, F0, stdDev, 1.0, delta);
    return bachelierBlackFormulaImpliedVol(Option::Call, K, F0, tte, price);
}

struct SkewFit {
    Real lambda = 0.0, delta = 0.0;
    bool saturated = false;
    Real residLow = 0.0, residHigh = 0.0; // target - model at the wing strikes, in normal vol
};

// match the anchor-strike vol exactly via lambda, and the (high - low) vol differential via delta
// within the admissible range [0, deltaMax]; on saturation, clamp and report residuals
SkewFit calibrateSkew(const Real F0, const Real Kanchor, const Real volAnchor, const Real Klow, const Real volLow,
                      const Real Khigh, const Real volHigh, const Real I, const Time tte, const Real deltaMax) {
    auto lambdaFor = [&](const Real delta) {
        const Real target = bachelierBlackFormula(Option::Call, Kanchor, F0, volAnchor * std::sqrt(tte));
        Brent b;
        return b.solve(
            [&](const Real lambda) {
                return blackFormula(Option::Call, Kanchor, F0, lambda * std::sqrt(I), 1.0, delta) - target;
            },
            1e-14, 0.01, 1e-8, 10.0);
    };
    auto skewFor = [&](const Real delta) {
        const Real lambda = lambdaFor(delta);
        return ddCapletNormalVol(F0, Khigh, lambda, delta, I, tte) -
               ddCapletNormalVol(F0, Klow, lambda, delta, I, tte);
    };
    const Real targetSkew = volHigh - volLow;
    const Real skewAtZero = skewFor(0.0), skewAtMax = skewFor(deltaMax);
    SkewFit fit;
    // skew is monotone decreasing in delta on [0, deltaMax] (asserted separately)
    if (targetSkew >= skewAtZero) {
        fit.delta = 0.0;
        fit.saturated = targetSkew > skewAtZero + 1e-12;
    } else if (targetSkew <= skewAtMax) {
        fit.delta = deltaMax;
        fit.saturated = targetSkew < skewAtMax - 1e-12;
    } else {
        Brent b;
        fit.delta = b.solve([&](const Real d) { return skewFor(d) - targetSkew; }, 1e-10, 0.5 * deltaMax, 0.0,
                            deltaMax);
    }
    fit.lambda = lambdaFor(fit.delta);
    fit.residLow = volLow - ddCapletNormalVol(F0, Klow, fit.lambda, fit.delta, I, tte);
    fit.residHigh = volHigh - ddCapletNormalVol(F0, Khigh, fit.lambda, fit.delta, I, tte);
    return fit;
}

} // namespace

BOOST_AUTO_TEST_CASE(testSkewShiftRoundTrip) {
    BOOST_TEST_MESSAGE("Skew-through-shift round trip: generate 3-strike vols from known (lambda*, delta*), "
                       "recover by ATM + skew calibration...");

    // unit-level parametrization: I_j from the exact integral engine with unit level
    FmmTestBed unitBed(FmmParametrization::LocalVolType::DisplacedDiffusion, 1.0, 0.02, 1, 0.0, 0.0);
    const Real F0 = 0.03;
    const std::vector<Size> buckets = {4, 12, 20};
    const std::vector<Real> trueDeltas = {0.0, 0.02, 2.0, 4.0}; // 4.0 = 1/tau (quarterly)
    for (const Size j : buckets) {
        const Real I = unitBed.parametrization->integratedCovariance(j, j, 0.0, unitBed.parametrization->rateTime(j));
        const Time tte = unitBed.parametrization->rateTime(j);
        for (const Real dStar : trueDeltas) {
            const Real lStar = 0.009 / (F0 + dStar); // ~90bp atm normal vol
            const Real Klow = F0 - 0.01, Khigh = F0 + 0.01;
            const Real vAtm = ddCapletNormalVol(F0, F0, lStar, dStar, I, tte);
            const Real vLow = ddCapletNormalVol(F0, Klow, lStar, dStar, I, tte);
            const Real vHigh = ddCapletNormalVol(F0, Khigh, lStar, dStar, I, tte);
            const SkewFit fit = calibrateSkew(F0, F0, vAtm, Klow, vLow, Khigh, vHigh, I, tte, 4.0);
            BOOST_CHECK_MESSAGE(!fit.saturated, "unexpected saturation for interior delta* " << dStar);
            // strong criterion: repricing of all three strikes to < 0.001 bp normal vol
            for (const Real K : {Klow, F0, Khigh}) {
                const Real vTgt = ddCapletNormalVol(F0, K, lStar, dStar, I, tte);
                const Real vFit = ddCapletNormalVol(F0, K, fit.lambda, fit.delta, I, tte);
                BOOST_CHECK_MESSAGE(std::fabs(vFit - vTgt) < 1e-7,
                                    "repricing failure j=" << j << " delta*=" << dStar << " K=" << K << ": "
                                                           << (vFit - vTgt) * 1e4 << " bp");
            }
            // parameter recovery where the skew is well-conditioned (away from the normal limit)
            if (dStar < 3.0)
                BOOST_CHECK_MESSAGE(std::fabs(fit.delta - dStar) < 1e-3,
                                    "delta recovery j=" << j << ": " << fit.delta << " vs " << dStar);
        }
    }
}

BOOST_AUTO_TEST_CASE(testSkewMonotoneAttainableRange) {
    BOOST_TEST_MESSAGE("Skew differential is monotone decreasing in delta (well-posed inversion); "
                       "documenting the attainable range...");
    FmmTestBed unitBed(FmmParametrization::LocalVolType::DisplacedDiffusion, 1.0, 0.02, 1, 0.0, 0.0);
    const Size j = 12;
    const Real F0 = 0.03, I = unitBed.parametrization->integratedCovariance(j, j, 0.0, 3.0);
    const Time tte = 3.0;
    Real prev = QL_MAX_REAL;
    for (Real d = 0.0; d <= 4.0 + 1e-12; d += 0.25) {
        const Real lambda = 0.009 / (F0 + d);
        const Real skew = ddCapletNormalVol(F0, 0.04, lambda, d, I, tte) -
                          ddCapletNormalVol(F0, 0.02, lambda, d, I, tte);
        BOOST_CHECK_MESSAGE(skew < prev + 1e-12, "skew not decreasing at delta=" << d);
        if (d == 0.0 || std::fabs(d - 4.0) < 1e-9)
            BOOST_TEST_MESSAGE("attainable skew endpoint at delta=" << d << ": " << skew * 1e4
                                                                    << " bp per 200bp strike width");
        prev = skew;
    }
}

BOOST_AUTO_TEST_CASE(testSkewSaturationOnMarketShapedData) {
    BOOST_TEST_MESSAGE("Skew calibration against real market-shaped USD-SOFR quotes: saturation and "
                       "curvature residuals must be reported, not absorbed...");
    // Quotes lifted from the pinned benchmark data (Examples/Products/Input/marketdata.csv,
    // asof 2025-02-10, CAPFLOOR/RATE_NVOL/USD/5Y/1D): flat cap normal vols used here as
    // market-shaped TARGETS for a single caplet bucket; proper cap->optionlet stripping is the
    // A3 calibration milestone's job.
    //   K=1%: 114.244 bp   K=3%: 105.458 bp   K=4%: 101.095 bp   K=5%: 110.116 bp
    FmmTestBed unitBed(FmmParametrization::LocalVolType::DisplacedDiffusion, 1.0, 0.02, 1, 0.0, 0.0);
    const Size j = 20; // 5y bucket
    const Real I = unitBed.parametrization->integratedCovariance(j, j, 0.0, 5.0);
    const Time tte = 5.0;
    const Real F0 = 0.042; // ~SOFR 5y forward area; anchor strike 4% is near-atm

    // (a) inverted (receiver) window K 1% -> 4%: target skew -13.1 bp is below the attainable
    //     minimum (~0 at delta = 1/tau) -> the calibrator must clamp at delta = 1/tau and report
    const SkewFit inv = calibrateSkew(F0, 0.04, 0.0101095, 0.01, 0.0114244, 0.04, 0.0101095, I, tte, 4.0);
    BOOST_TEST_MESSAGE("(a) inverted window: delta=" << inv.delta << " saturated=" << inv.saturated
                                                     << " residual at K=1%: " << inv.residLow * 1e4 << " bp");
    BOOST_CHECK_MESSAGE(inv.saturated && std::fabs(inv.delta - 4.0) < 1e-12,
                        "expected saturation at delta = 1/tau, got delta=" << inv.delta);
    BOOST_CHECK_MESSAGE(inv.residLow > 0.0005,
                        "expected a reported unattainable-skew residual > 5bp, got " << inv.residLow * 1e4 << " bp");
    // anchor strike must still reprice exactly
    BOOST_CHECK_SMALL(ddCapletNormalVol(F0, 0.04, inv.lambda, inv.delta, I, tte) - 0.0101095, 1e-10);

    // (b) upper window K 3% -> 5% around the smile minimum: the +4.66 bp slope is attainable, but
    //     the V-shape (curvature) is not a DD degree of freedom -> wings reprice, the smile bottom
    //     at 4% shows the reported curvature residual
    const SkewFit vfit = calibrateSkew(F0, 0.04, 0.0101095, 0.03, 0.0105458, 0.05, 0.0110116, I, tte, 4.0);
    const Real modelAt3 = ddCapletNormalVol(F0, 0.03, vfit.lambda, vfit.delta, I, tte);
    const Real modelAt5 = ddCapletNormalVol(F0, 0.05, vfit.lambda, vfit.delta, I, tte);
    BOOST_TEST_MESSAGE("(b) smile window: delta=" << vfit.delta << " saturated=" << vfit.saturated
                                                  << "; residuals: K=3% " << (0.0105458 - modelAt3) * 1e4
                                                  << " bp, K=5% " << (0.0110116 - modelAt5) * 1e4 << " bp");
    // note: anchoring the smile BOTTOM while matching the wing differential leaves symmetric
    // wing residuals of equal sign - the honest signature of missing curvature
    BOOST_CHECK_MESSAGE(!vfit.saturated, "slope within attainable family should not saturate");
    BOOST_CHECK_SMALL((0.0105458 - modelAt3) - (0.0110116 - modelAt5), 5e-5); // near-equal wing residuals
    BOOST_CHECK_MESSAGE(0.0105458 - modelAt3 > 0.0002,
                        "expected reported curvature residual > 2bp at the wings, got "
                            << (0.0105458 - modelAt3) * 1e4 << " bp");
}

BOOST_AUTO_TEST_CASE(testSwaptionApproxAccuracyTable) {
    BOOST_TEST_MESSAGE("A3 acceptance 1: frozen-gradient swaption approximations (strike-flat normal and "
                       "effective-shift) vs same-model MC over an expiry x tenor x strike grid; randomized-"
                       "Sobol replications for ~0.1bp reference uncertainty...");

    struct Cfg {
        FmmParametrization::LocalVolType type;
        Real level, shift;
        const char* label;
    };
    const std::vector<Cfg> cfgs = {{FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, "DD 1/tau"},
                                   {FmmParametrization::LocalVolType::DisplacedDiffusion, 0.15, 0.02, "DD 2%"}};

    const std::vector<Size> expiries = {4, 8, 12}; // 1y, 2y, 3y
    const std::vector<Size> tenors = {4, 8};       // 1y, 2y
    const std::vector<Real> offsets = {-0.01, -0.005, 0.0, 0.005, 0.01};
    const std::vector<const char*> kLab = {"-100", "-50", "atm", "+50", "+100"};
    const Size reps = 8, pathsPerRep = 16384;

    struct Worst {
        Real bp = 0.0;
        std::string cell;
        void track(const Real d, const std::string& c) {
            if (std::fabs(d) > std::fabs(bp)) {
                bp = d;
                cell = c;
            }
        }
    } worstAtmFlat, worstAtmShift, worstWingFlat, worstWingShift;

    for (const auto& cfg : cfgs) {
        FmmTestBed bed(cfg.type, cfg.level, cfg.shift, 3);
        const auto& p = *bed.parametrization;
        for (const Size a : expiries) {
            const Time Te = p.rateTime(a);
            // one simulation batch per expiry serves every (tenor, strike) cell
            std::vector<FmmSwapSpec> swaps;
            std::vector<Real> strikes;
            for (const Size len : tenors) {
                FmmSwapSpec s;
                s.a = a;
                s.b = a + len;
                for (Size c = a + 4; c <= s.b; c += 4) {
                    s.fixedPayIndices.push_back(c);
                    s.fixedAccruals.push_back(1.0);
                }
                const Real S0 = fmmForwardSwapRate(p, s);
                for (const Real off : offsets) {
                    swaps.push_back(s);
                    strikes.push_back(S0 + off);
                }
            }
            const auto mc = fmmSwaptionMc(bed.model, swaps, strikes, Option::Call, pathsPerRep, reps,
                                          7000 + 100 * static_cast<BigNatural>(a));
            for (Size f = 0; f < swaps.size(); ++f) {
                const Size ki = f % offsets.size();
                const auto flat = fmmSwaptionApprox(p, swaps[f], strikes[f]);
                const auto shift =
                    fmmSwaptionApprox(p, swaps[f], strikes[f], Option::Call, FmmSwaptionApproxMethod::EffectiveShift);
                const Real seBp = mc[f].normalVolSe * 1e4;
                const Real dFlat = (flat.normalVol - mc[f].normalVol) * 1e4;
                const Real dShift = (shift.normalVol - mc[f].normalVol) * 1e4;
                std::ostringstream cell;
                cell << cfg.label << " " << Te << "y x " << 0.25 * (swaps[f].b - swaps[f].a) << "y " << kLab[ki];
                BOOST_TEST_MESSAGE("APPROX ROW | " << cell.str() << " | mc " << mc[f].normalVol * 1e4 << " +/- " << seBp
                                                   << " bp | strike-flat " << flat.normalVol * 1e4 << " (diff " << dFlat
                                                   << ") | effective-shift " << shift.normalVol * 1e4 << " (diff "
                                                   << dShift << ", delta_S " << shift.effectiveShift << ")");
                if (ki == 2) {
                    worstAtmFlat.track(dFlat, cell.str());
                    worstAtmShift.track(dShift, cell.str());
                    BOOST_CHECK_MESSAGE(std::fabs(dFlat) < 1.0 + 3.0 * seBp,
                                        cell.str() << ": strike-flat ATM error " << dFlat << " bp");
                    BOOST_CHECK_MESSAGE(std::fabs(dShift) < 1.0 + 3.0 * seBp,
                                        cell.str() << ": effective-shift ATM error " << dShift << " bp");
                } else {
                    worstWingFlat.track(dFlat, cell.str());
                    worstWingShift.track(dShift, cell.str());
                    // strike-flat wings are reported (known approximation error); the effective
                    // shift carries the deal-strike accuracy claim
                    BOOST_CHECK_MESSAGE(std::fabs(dShift) < 2.0 + 3.0 * seBp,
                                        cell.str() << ": effective-shift wing error " << dShift << " bp");
                }
            }
        }
    }
    BOOST_TEST_MESSAGE("WORST ATM: strike-flat " << worstAtmFlat.bp << " bp (" << worstAtmFlat.cell << "), effective-shift "
                                                 << worstAtmShift.bp << " bp (" << worstAtmShift.cell << ")");
    BOOST_TEST_MESSAGE("WORST WING: strike-flat " << worstWingFlat.bp << " bp (" << worstWingFlat.cell
                                                  << "), effective-shift " << worstWingShift.bp << " bp ("
                                                  << worstWingShift.cell << ")");
}

namespace {

// separable testbed: M=20 quarterly grid, DD shift (default 1/tau), a(t) segments at {1y, 2y, 3y}
struct CalibBed {
    explicit CalibBed(const Real flatRate, const Real rhoInf = 0.6, const Real shift = 4.0, const Real level = 0.002) {
        Settings::instance().evaluationDate() = Date(19, September, 2026);
        curve = Handle<YieldTermStructure>(
            QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), flatRate, Actual365Fixed()));
        const Size M = 20;
        Array rateTimes(M + 1);
        for (Size k = 0; k <= M; ++k)
            rateTimes[k] = 0.25 * static_cast<Real>(k);
        Array shifts(M, shift);
        segTimes = Array(3);
        segTimes[0] = 1.0;
        segTimes[1] = 2.0;
        segTimes[2] = 3.0;
        std::vector<Array> volLevels(M, Array(4, level));
        parametrization = QuantLib::ext::make_shared<FmmParametrization>(
            EURCurrency(), curve, rateTimes, shifts, segTimes, volLevels,
            FmmParametrization::LocalVolType::DisplacedDiffusion, rhoInf, 0.08, 3);
    }
    FmmSwapSpec coterminal(const Size a, const Size b = 20) const {
        FmmSwapSpec s;
        s.a = a;
        s.b = b;
        for (Size c = a + 4; c <= b; c += 4) {
            s.fixedPayIndices.push_back(c);
            s.fixedAccruals.push_back(1.0);
        }
        return s;
    }
    Real atmCapletVol(const Size j) const {
        const Real F0 = (curve->discount(parametrization->rateTime(j - 1)) /
                             curve->discount(parametrization->rateTime(j)) -
                         1.0) /
                        parametrization->tau(j);
        return fmmCapletNormalVol(*parametrization, j, F0, true);
    }
    Handle<YieldTermStructure> curve;
    Array segTimes;
    QuantLib::ext::shared_ptr<FmmParametrization> parametrization;
};

} // namespace

BOOST_AUTO_TEST_CASE(testCalibrationRoundTrip) {
    BOOST_TEST_MESSAGE("A3 acceptance 2: joint bootstrap round trip - generate targets from known "
                       "(Lambda*, a*), recover them...");

    CalibBed bed(0.03);
    // true separable state
    FmmSeparableVols truth;
    truth.segmentTimes = bed.segTimes;
    truth.a = {1.0, 1.15, 0.9, 1.05};
    for (Size j = 1; j <= 20; ++j)
        truth.levels.push_back(0.0022 + 0.00002 * static_cast<Real>(j));
    truth.apply(*bed.parametrization);

    // targets generated from the model
    std::vector<FmmCapletVolTarget> capletTargets;
    for (Size j = 1; j <= 20; ++j)
        capletTargets.push_back({j, bed.atmCapletVol(j), true});
    std::vector<FmmSwaptionVolTarget> swaptionTargets;
    for (const Size a : {4, 8, 12, 16}) {
        FmmSwaptionVolTarget t;
        t.swap = bed.coterminal(a);
        t.normalVol = fmmSwaptionApprox(*bed.parametrization, t.swap,
                                        fmmForwardSwapRate(*bed.parametrization, t.swap))
                          .normalVol;
        std::ostringstream lbl;
        lbl << "cot_" << a / 4 << "y";
        t.label = lbl.str();
        swaptionTargets.push_back(t);
    }

    // calibrate from a cold start
    FmmSeparableVols v;
    v.segmentTimes = bed.segTimes;
    v.a = {1.0, 1.0, 1.0, 1.0};
    v.levels.assign(20, 0.004);
    const auto report = fmmJointBootstrap(*bed.parametrization, v, capletTargets, swaptionTargets);
    BOOST_TEST_MESSAGE("round trip: iterations " << report.iterations << ", runtime " << report.runtimeSeconds
                                                 << " s, max |error| " << report.maxAbsErrorBp() << " bp");
    BOOST_CHECK_MESSAGE(report.converged, "joint bootstrap did not converge");
    BOOST_CHECK_MESSAGE(report.maxAbsErrorBp() < 1e-3, "repricing error " << report.maxAbsErrorBp() << " bp");
    // parameter recovery: the effective lambda_j(t) products must match the truth
    Real worstRel = 0.0;
    for (Size j = 1; j <= 20; ++j)
        for (Size k = 0; k < 4; ++k) {
            const Time tMid = (k == 0 ? 0.5 : bed.segTimes[k - 1] + 0.5);
            if (tMid >= bed.parametrization->rateTime(j))
                continue;
            const Real fit = bed.parametrization->volLevel(j, tMid);
            const Real tru = truth.levels[j - 1] * truth.a[k];
            worstRel = std::max(worstRel, std::fabs(fit / tru - 1.0));
        }
    BOOST_TEST_MESSAGE("worst relative parameter error " << worstRel);
    BOOST_CHECK_MESSAGE(worstRel < 5e-5, "parameter recovery error " << worstRel); // solver acc 1e-8 in vol space
    BOOST_CHECK_SMALL(v.a[0] - 1.0, 1e-14); // identifiability normalization
}

BOOST_AUTO_TEST_CASE(testCalibrationOnMarketQuotes) {
    BOOST_TEST_MESSAGE("A3 acceptance 3 (QuantExt level): joint calibration to real USD-SOFR quotes "
                       "(Products marketdata.csv, 2025-02-10) - convergence, residuals, runtime, "
                       "bump stability...");
    // Provenance: SWAPTION/RATE_NVOL/USD/SOFR coterminal-to-5y ATM quotes and CAPFLOOR/RATE_NVOL
    // USD ATM cap vols used as per-bucket caplet-vol proxies (flat-cap-to-optionlet stripping is
    // the A5/ORE-plumbing job; the mapping choice is documented in the calibration report).
    // Curve: flat 4.20% stands in at the observed SOFR level; curve-accurate bootstrap in A5.
    auto runCalibration = [](const Real flatRate, const Real rhoInf, FmmSeparableVols& v,
                             FmmCalibrationReport& report) {
        CalibBed bed(flatRate, rhoInf);
        std::vector<FmmCapletVolTarget> capletTargets = {{4, 0.0055467, true},
                                                         {8, 0.0083185, true},
                                                         {12, 0.0094564, true},
                                                         {16, 0.0099147, true},
                                                         {20, 0.0101515, true}};
        std::vector<FmmSwaptionVolTarget> swaptionTargets;
        const std::vector<std::pair<Size, Real>> cot = {
            {4, 0.01029015}, {8, 0.0105039}, {12, 0.0105350}, {16, 0.0104501}};
        for (const auto& c : cot) {
            FmmSwaptionVolTarget t;
            t.swap = bed.coterminal(c.first);
            t.normalVol = c.second;
            std::ostringstream lbl;
            lbl << "cot_" << c.first / 4 << "y_to_5y";
            t.label = lbl.str();
            swaptionTargets.push_back(t);
        }
        v.segmentTimes = bed.segTimes;
        v.a = {1.0, 1.0, 1.0, 1.0};
        v.levels.assign(20, 0.0025);
        report = fmmJointBootstrap(*bed.parametrization, v, capletTargets, swaptionTargets, 15, 0.01);
    };

    // ---- part 1: strategy (b) on real coterminal swaption quotes only - exact fit expected ----
    auto runSwaptionOnly = [](const Real flatRate, FmmSeparableVols& v) {
        CalibBed bed(flatRate, 0.9);
        std::vector<FmmSwaptionVolTarget> swaptionTargets;
        const std::vector<std::pair<Size, Real>> cot = {
            {4, 0.01029015}, {8, 0.0105039}, {12, 0.0105350}, {16, 0.0104501}};
        for (const auto& c : cot) {
            FmmSwaptionVolTarget t;
            t.swap = bed.coterminal(c.first);
            t.normalVol = c.second;
            swaptionTargets.push_back(t);
        }
        v.segmentTimes = bed.segTimes;
        v.a = {1.0, 1.0, 1.0, 1.0};
        v.levels.assign(20, 0.0025);
        fmmSwaptionTimeDependenceBootstrap(*bed.parametrization, v, swaptionTargets);
        Real worst = 0.0;
        for (const auto& t : swaptionTargets)
            worst = std::max(worst, std::fabs(fmmSwaptionApprox(*bed.parametrization, t.swap,
                                                                fmmForwardSwapRate(*bed.parametrization, t.swap))
                                                  .normalVol -
                                              t.normalVol) *
                                        1e4);
        return worst;
    };
    FmmSeparableVols vB;
    const Real worstB = runSwaptionOnly(0.042, vB);
    BOOST_TEST_MESSAGE("strategy (b), real coterminals: worst residual " << worstB << " bp; a(t) = {"
                                                                         << vB.a[0] << ", " << vB.a[1] << ", "
                                                                         << vB.a[2] << ", " << vB.a[3] << "}");
    BOOST_CHECK_MESSAGE(worstB < 0.01, "strategy (b) residual " << worstB << " bp");
    BOOST_CHECK_SMALL(vB.a[0] - 1.0, 1e-14);

    // ---- part 2: joint (c) with FLAT-CAP vols as caplet proxies - a PROXY-ARTIFACT record
    // (owner correction 2026-09-19): the residual below reflects the proxy's inconsistency (flat
    // cap vols are not optionlet vols; they average the small front optionlets of the steep 2025
    // vol term structure), NOT an economic caplet/swaption basis. Correlation cannot close it
    // (mean residual flat in rhoInf). Real-market JOINT calibration is PENDING until validated
    // against correctly stripped SOFR optionlets or actual cap prices (A5 stripping, fail-fast).
    FmmSeparableVols v;
    FmmCalibrationReport report;
    runCalibration(0.042, 0.9, v, report);
    Real worstCapletBp = 0.0, worstSwaptionBp = 0.0;
    for (const auto& r : report.rows) {
        BOOST_TEST_MESSAGE(r.instrument << ": market " << r.marketVol * 1e4 << " bp, model " << r.modelVol * 1e4
                                        << " bp, error " << r.errorBp << " bp");
        if (r.instrument.rfind("caplet", 0) == 0)
            worstCapletBp = std::max(worstCapletBp, std::fabs(r.errorBp));
        else
            worstSwaptionBp = std::max(worstSwaptionBp, std::fabs(r.errorBp));
    }
    BOOST_TEST_MESSAGE("joint on cap-proxy targets: iterations " << report.iterations << ", runtime "
                                                                 << report.runtimeSeconds << " s; swaptions exact to "
                                                                 << worstSwaptionBp << " bp, cap-proxy basis up to "
                                                                 << worstCapletBp << " bp");
    BOOST_CHECK_MESSAGE(worstSwaptionBp < 0.01, "swaption residual " << worstSwaptionBp << " bp");
    BOOST_CHECK_MESSAGE(worstCapletBp < 5.0, "cap-proxy basis unexpectedly large: " << worstCapletBp << " bp");

    // ---- part 3: correlation recovery on CONSISTENT targets - generate caplet + swaption vols
    // from rhoInf* = 0.85, recalibrate with the outer solve, recover the correlation ----
    auto makeConsistentTargets = [&](const Real rhoInf, std::vector<FmmCapletVolTarget>& caps,
                                     std::vector<FmmSwaptionVolTarget>& swps) {
        CalibBed bed(0.042, rhoInf);
        FmmSeparableVols truth;
        truth.segmentTimes = bed.segTimes;
        truth.a = {1.0, 1.1, 0.95, 1.02};
        truth.levels.assign(20, 0.0026);
        truth.apply(*bed.parametrization);
        for (const Size j : {4, 8, 12, 16, 20})
            caps.push_back({j, bed.atmCapletVol(j), true});
        for (const Size a : {4, 8, 12, 16}) {
            FmmSwaptionVolTarget t;
            t.swap = bed.coterminal(a);
            t.normalVol = fmmSwaptionApprox(*bed.parametrization, t.swap,
                                            fmmForwardSwapRate(*bed.parametrization, t.swap))
                              .normalVol;
            swps.push_back(t);
        }
    };
    std::vector<FmmCapletVolTarget> capsC;
    std::vector<FmmSwaptionVolTarget> swpsC;
    makeConsistentTargets(0.85, capsC, swpsC);
    auto meanCapletResidual = [&](const Real rhoInf) {
        CalibBed bed(0.042, rhoInf);
        FmmSeparableVols vT;
        vT.segmentTimes = bed.segTimes;
        vT.a = {1.0, 1.0, 1.0, 1.0};
        vT.levels.assign(20, 0.0025);
        const auto rep = fmmJointBootstrap(*bed.parametrization, vT, capsC, swpsC, 12, 0.01);
        Real sum = 0.0;
        Size n = 0;
        for (const auto& r : rep.rows)
            if (r.instrument.rfind("caplet", 0) == 0) {
                sum += r.errorBp;
                ++n;
            }
        return sum / static_cast<Real>(n);
    };
    Brent rhoSolver;
    const Real rhoStar = rhoSolver.solve(meanCapletResidual, 1e-4, 0.8, 0.40, 0.999);
    BOOST_TEST_MESSAGE("correlation recovery: rhoInf* = " << rhoStar << " (truth 0.85)");
    BOOST_CHECK_MESSAGE(std::fabs(rhoStar - 0.85) < 0.02, "rho recovery failed: " << rhoStar);

    // ---- owner requirement: parameter stability under a +10 bp parallel bump (strategy (b)) ----
    FmmSeparableVols vBumped;
    const Real worstBumped = runSwaptionOnly(0.043, vBumped);
    Real worstMove = 0.0;
    for (Size j = 0; j < 20; ++j)
        for (Size k = 0; k < 4; ++k)
            worstMove = std::max(worstMove, std::fabs(vBumped.levels[j] * vBumped.a[k] /
                                                          (vB.levels[j] * vB.a[k]) -
                                                      1.0));
    BOOST_TEST_MESSAGE("bump stability: worst relative parameter move " << worstMove << ", bumped residual "
                                                                        << worstBumped << " bp");
    BOOST_CHECK_MESSAGE(worstMove < 0.10, "unstable under bump: " << worstMove);
    BOOST_CHECK_MESSAGE(worstBumped < 0.01, "bumped fit degraded: " << worstBumped << " bp");
}

BOOST_AUTO_TEST_CASE(testReplicationRepricesLgmSwaptions) {
    BOOST_TEST_MESSAGE("A3 acceptance 4: replication-mode FMM fed LGM's parameters reprices LGM's "
                       "calibration swaptions (vs AnalyticLgmSwaptionEngine), schedule-accurate...");

    const Date asof(19, September, 2026);
    Settings::instance().evaluationDate() = asof;
    Handle<YieldTermStructure> curve(
        QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.03, Actual365Fixed()));
    const Actual365Fixed dc;

    // quarterly date grid to 5y, NullCalendar/Unadjusted so schedule and model times align exactly
    const Date start = asof;
    const Date end = asof + Period(5, Years);
    const Schedule quarterly(start, end, Period(3, Months), NullCalendar(), Unadjusted, Unadjusted,
                             DateGeneration::Forward, false);
    const Size M = quarterly.size() - 1;
    Array rateTimes(M + 1);
    for (Size k = 0; k <= M; ++k)
        rateTimes[k] = dc.yearFraction(asof, quarterly[k]);

    // LGM with structured piecewise-constant alpha (breaks at 1y, 2y, 3y) and constant kappa
    Array alphaTimes(3);
    alphaTimes[0] = rateTimes[4];
    alphaTimes[1] = rateTimes[8];
    alphaTimes[2] = rateTimes[12];
    Array alpha(4);
    alpha[0] = 0.0090;
    alpha[1] = 0.0110;
    alpha[2] = 0.0100;
    alpha[3] = 0.0095;
    auto lgmParam = QuantLib::ext::make_shared<IrLgm1fPiecewiseConstantParametrization>(
        EURCurrency(), curve, alphaTimes, alpha, Array(), Array(1, 0.01));
    auto lgmModel = QuantLib::ext::make_shared<LinearGaussMarkovModel>(lgmParam);

    // replication-mode FMM on the schedule grid, fed the SAME parametrization
    auto fmmParam = QuantLib::ext::make_shared<FmmParametrization>(EURCurrency(), curve, rateTimes, lgmParam);

    auto index = QuantLib::ext::make_shared<IborIndex>("FMMTEST", Period(3, Months), 0, EURCurrency(),
                                                       NullCalendar(), Unadjusted, false, dc, curve);
    auto swapEngine = QuantLib::ext::make_shared<DiscountingSwapEngine>(curve);
    auto lgmEngine = QuantLib::ext::make_shared<AnalyticLgmSwaptionEngine>(lgmParam, curve);

    Real worstBp = 0.0;
    for (const Size a : {4, 8, 12, 16}) {
        const Date expiry = quarterly[a];
        const Schedule fixedSched(expiry, end, Period(1, Years), NullCalendar(), Unadjusted, Unadjusted,
                                  DateGeneration::Forward, false);
        const Schedule floatSched(expiry, end, Period(3, Months), NullCalendar(), Unadjusted, Unadjusted,
                                  DateGeneration::Forward, false);
        // ATM strike from the discounting engine
        VanillaSwap probe(VanillaSwap::Payer, 1.0, fixedSched, 0.03, dc, floatSched, index, 0.0, dc);
        probe.setPricingEngine(swapEngine);
        const Real fair = probe.fairRate();
        auto underlying = QuantLib::ext::make_shared<VanillaSwap>(VanillaSwap::Payer, 1.0, fixedSched, fair, dc,
                                                                  floatSched, index, 0.0, dc);
        underlying->setPricingEngine(swapEngine);
        auto swaption = QuantLib::ext::make_shared<Swaption>(
            underlying, QuantLib::ext::make_shared<EuropeanExercise>(expiry));
        swaption->setPricingEngine(lgmEngine);
        const Real lgmPrice = swaption->NPV();

        // grid-mapped FMM swap spec with the schedule's actual accruals on the fixed leg
        FmmSwapSpec spec;
        spec.a = a;
        spec.b = M;
        for (Size c = a + 4; c <= M; c += 4) {
            spec.fixedPayIndices.push_back(c);
            spec.fixedAccruals.push_back(rateTimes[c] - rateTimes[c - 4]);
        }
        const Real S0 = fmmForwardSwapRate(*fmmParam, spec);
        BOOST_REQUIRE_MESSAGE(std::fabs(S0 - fair) < 1e-10,
                              "forward swap rate mismatch: fmm " << S0 << " vs discounting " << fair);
        const auto apx = fmmSwaptionApprox(*fmmParam, spec, fair);

        const Time Te = rateTimes[a];
        const Real A0 = fmmAnnuity(*fmmParam, spec);
        const Real lgmVol = bachelierBlackFormulaImpliedVol(Option::Call, fair, S0, Te, lgmPrice / A0);
        const Real diffBp = (apx.normalVol - lgmVol) * 1e4;
        BOOST_TEST_MESSAGE("coterminal " << a / 4 << "y x " << (M - a) / 4 << "y: LGM " << lgmVol * 1e4
                                         << " bp, FMM-replication approx " << apx.normalVol * 1e4 << " bp, diff "
                                         << diffBp << " bp");
        worstBp = std::max(worstBp, std::fabs(diffBp));
        BOOST_CHECK_MESSAGE(std::fabs(diffBp) < 1.0,
                            "replication repricing error " << diffBp << " bp at expiry index " << a);
    }
    BOOST_TEST_MESSAGE("worst replication repricing error: " << worstBp << " bp (approximation-level target < 1 bp)");
}

namespace {

// payer swap flows into an FmmCallableInstrument flow set: receive float quarterly, pay fixed K
// annually (grid indices multiples of 4), unit notional, from the FLOAT RECEIVER's perspective
void fillPayerSwapFlows(FmmCallableInstrument& inst, const FmmParametrization& p, const Size from, const Size to,
                        const Real K, const Real sign = 1.0) {
    inst.lastFlowIdx = to;
    inst.fixedFlows.assign(p.numberOfRates() + 1, 0.0);
    inst.floatWeights.assign(p.numberOfRates() + 1, 0.0);
    for (Size j = from + 1; j <= to; ++j)
        inst.floatWeights[j] = sign;
    for (Size c = from + 4; c <= to; c += 4)
        inst.fixedFlows[c] = -sign * K * (p.rateTime(c) - p.rateTime(c - 4));
}

// replication-mode bed on real dates: quarterly grid to 5y on a flat 3% curve, piecewise-constant
// LGM alpha (the setup of testLsmBermudanVsLgmGrid), FMM built as the exact replication of the LGM
struct ReplicationBed {
    ReplicationBed()
        : asof(19, September, 2026), end(asof + Period(5, Years)),
          quarterly(asof, end, Period(3, Months), NullCalendar(), Unadjusted, Unadjusted, DateGeneration::Forward,
                    false) {
        Settings::instance().evaluationDate() = asof;
        curve = Handle<YieldTermStructure>(
            QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.03, Actual365Fixed()));
        M = quarterly.size() - 1;
        Array rateTimes(M + 1);
        for (Size k = 0; k <= M; ++k)
            rateTimes[k] = dc.yearFraction(asof, quarterly[k]);
        Array alphaTimes(3);
        alphaTimes[0] = rateTimes[4];
        alphaTimes[1] = rateTimes[8];
        alphaTimes[2] = rateTimes[12];
        Array alpha(4);
        alpha[0] = 0.0090;
        alpha[1] = 0.0110;
        alpha[2] = 0.0100;
        alpha[3] = 0.0095;
        lgmParam = QuantLib::ext::make_shared<IrLgm1fPiecewiseConstantParametrization>(EURCurrency(), curve, alphaTimes,
                                                                                       alpha, Array(), Array(1, 0.01));
        lgmModel = QuantLib::ext::make_shared<LinearGaussMarkovModel>(lgmParam);
        fmmParam = QuantLib::ext::make_shared<FmmParametrization>(EURCurrency(), curve, rateTimes, lgmParam);
        fmmModel = QuantLib::ext::make_shared<ForwardMarketModel>(fmmParam);
    }
    Date asof, end;
    Schedule quarterly;
    Actual365Fixed dc;
    Handle<YieldTermStructure> curve;
    Size M = 0;
    QuantLib::ext::shared_ptr<IrLgm1fPiecewiseConstantParametrization> lgmParam;
    QuantLib::ext::shared_ptr<LinearGaussMarkovModel> lgmModel;
    QuantLib::ext::shared_ptr<FmmParametrization> fmmParam;
    QuantLib::ext::shared_ptr<ForwardMarketModel> fmmModel;
};

// paired common-random-number identity Cancel == underlying + Enter: both pricers are valued on
// the SAME paths (same seed and sequence), the realized underlying total is replaced path by path
// by its curve value, and the pathwise difference is reported with its confidence interval. The
// FMM's bank-account / floating-flow identities make the difference vanish on every path where
// the two policies decide alike, so this isolates policy-decision differences. Absolute
// tolerance per the revised validation contract, with the precision required to resolve it.
void pairedIdentityCheck(const FmmLsmPricer& cancelPricer, const FmmLsmPricer& enterPricer,
                         const Real underlyingCurve, const std::string& label, const Real optionValue,
                         const BigNatural seed = 31337) {
    std::vector<Real> vc, tc, ve, te;
    cancelPricer.pathValues(seed, vc, tc);
    enterPricer.pathValues(seed, ve, te);
    BOOST_REQUIRE_EQUAL(vc.size(), ve.size());
    IncrementalStatistics d;
    for (Size n = 0; n < vc.size(); ++n)
        d.add((vc[n] - tc[n] + underlyingCurve) - (underlyingCurve + ve[n]));
    // absolute tolerance, unit notional: 0.5% of the embedded option value (floor 2e-5)
    const Real tol = std::max(0.005 * std::fabs(optionValue), 2e-5);
    BOOST_TEST_MESSAGE("PAIRED-CRN IDENTITY | " << label << " | Cancel - underlying - Enter on common paths: "
                                                << d.mean() << " +/- " << d.errorEstimate() << " (95% CI ["
                                                << d.mean() - 1.96 * d.errorEstimate() << ", "
                                                << d.mean() + 1.96 * d.errorEstimate() << "]), absolute tolerance "
                                                << tol << " (0.5% of the option value " << optionValue << ")");
    BOOST_CHECK_MESSAGE(d.errorEstimate() < tol / 3.0,
                        label << ": paired precision " << d.errorEstimate() << " cannot resolve the tolerance " << tol);
    BOOST_CHECK_MESSAGE(std::fabs(d.mean()) < tol, label << ": paired identity difference " << d.mean() << " exceeds " << tol);
}

} // namespace

BOOST_AUTO_TEST_CASE(testSwaptionApproxWingsVsLgmAnalytic) {
    BOOST_TEST_MESSAGE("A3 strike-aware approximation, noise-free reference: replication-mode FMM vs "
                       "AnalyticLgmSwaptionEngine at ATM and +/-50, +/-100 bp strikes (coterminals 1y..4y "
                       "into 5y), strike-flat and effective-shift...");
    ReplicationBed rb;
    const auto& p = *rb.fmmParam;
    auto index = QuantLib::ext::make_shared<IborIndex>("FMMTEST", Period(3, Months), 0, EURCurrency(), NullCalendar(),
                                                       Unadjusted, false, rb.dc, rb.curve);
    auto swapEngine = QuantLib::ext::make_shared<DiscountingSwapEngine>(rb.curve);
    auto lgmEngine = QuantLib::ext::make_shared<AnalyticLgmSwaptionEngine>(rb.lgmParam, rb.curve);
    const std::vector<Real> offsets = {-0.01, -0.005, 0.0, 0.005, 0.01};
    const std::vector<const char*> kLab = {"-100", "-50", "atm", "+50", "+100"};
    Real worstAtmFlat = 0.0, worstAtmShift = 0.0, worstWingFlat = 0.0, worstWingShift = 0.0;
    for (const Size a : {4, 8, 12, 16}) {
        const Date expiry = rb.quarterly[a];
        const Schedule fixedSched(expiry, rb.end, Period(1, Years), NullCalendar(), Unadjusted, Unadjusted,
                                  DateGeneration::Forward, false);
        const Schedule floatSched(expiry, rb.end, Period(3, Months), NullCalendar(), Unadjusted, Unadjusted,
                                  DateGeneration::Forward, false);
        VanillaSwap probe(VanillaSwap::Payer, 1.0, fixedSched, 0.03, rb.dc, floatSched, index, 0.0, rb.dc);
        probe.setPricingEngine(swapEngine);
        const Real fair = probe.fairRate();
        FmmSwapSpec spec;
        spec.a = a;
        spec.b = rb.M;
        for (Size c = a + 4; c <= rb.M; c += 4) {
            spec.fixedPayIndices.push_back(c);
            spec.fixedAccruals.push_back(p.rateTime(c) - p.rateTime(c - 4));
        }
        const Real S0 = fmmForwardSwapRate(p, spec);
        BOOST_REQUIRE_MESSAGE(std::fabs(S0 - fair) < 1e-10, "forward swap rate mismatch: " << S0 - fair);
        const Time Te = p.rateTime(a);
        const Real A0 = fmmAnnuity(p, spec);
        for (Size ki = 0; ki < offsets.size(); ++ki) {
            const Real K = fair + offsets[ki];
            auto underlying = QuantLib::ext::make_shared<VanillaSwap>(VanillaSwap::Payer, 1.0, fixedSched, K, rb.dc,
                                                                      floatSched, index, 0.0, rb.dc);
            underlying->setPricingEngine(swapEngine);
            auto swaption =
                QuantLib::ext::make_shared<Swaption>(underlying, QuantLib::ext::make_shared<EuropeanExercise>(expiry));
            swaption->setPricingEngine(lgmEngine);
            const Real lgmVol = bachelierBlackFormulaImpliedVol(Option::Call, K, S0, Te, swaption->NPV() / A0);
            const auto flat = fmmSwaptionApprox(p, spec, K);
            const auto shift = fmmSwaptionApprox(p, spec, K, Option::Call, FmmSwaptionApproxMethod::EffectiveShift);
            const Real dFlat = (flat.normalVol - lgmVol) * 1e4, dShift = (shift.normalVol - lgmVol) * 1e4;
            BOOST_TEST_MESSAGE("LGM-ANALYTIC ROW | coterminal " << a / 4 << "y x " << (rb.M - a) / 4 << "y " << kLab[ki]
                                                                << " | LGM " << lgmVol * 1e4 << " bp | strike-flat diff "
                                                                << dFlat << " | effective-shift diff " << dShift
                                                                << " (delta_S " << shift.effectiveShift << ")");
            if (ki == 2) {
                worstAtmFlat = std::max(worstAtmFlat, std::fabs(dFlat));
                worstAtmShift = std::max(worstAtmShift, std::fabs(dShift));
            } else {
                worstWingFlat = std::max(worstWingFlat, std::fabs(dFlat));
                worstWingShift = std::max(worstWingShift, std::fabs(dShift));
            }
            BOOST_CHECK_MESSAGE(std::fabs(dShift) < 1.0, "effective-shift vs LGM analytic: " << dShift << " bp at "
                                                                                             << kLab[ki]);
        }
    }
    BOOST_TEST_MESSAGE("worst vs LGM analytic: ATM strike-flat " << worstAtmFlat << " bp, effective-shift "
                                                                 << worstAtmShift << " bp; wings strike-flat "
                                                                 << worstWingFlat << " bp, effective-shift "
                                                                 << worstWingShift << " bp");
}

BOOST_AUTO_TEST_CASE(testDealStrikeCalibration) {
    BOOST_TEST_MESSAGE("A3 strike-aware calibration: strategy (b) at DEAL strikes (30 bp below ATM, a Hagan-"
                       "style equivalent-spread offset) with the effective-shift approximation on the real "
                       "USD-SOFR coterminal quotes; MC-corrected loop (randomized-Sobol replications) under "
                       "DD 2% dynamics...");
    // the pinned market set carries ATM quotes only: the ATM normal vol is used at the deal strike
    // (documented simplification - no smile data), which exercises the strike-aware machinery
    // without asserting a market smile
    const std::vector<std::pair<Size, Real>> cot = {{4, 0.01029015}, {8, 0.0105039}, {12, 0.0105350}, {16, 0.0104501}};
    struct Cfg {
        Real shift, level;
        const char* label;
        bool mcLoop;
    };
    const std::vector<Cfg> cfgs = {{4.0, 0.002, "DD 1/tau", false}, {0.02, 0.15, "DD 2%", true}};
    for (const auto& cfg : cfgs) {
        CalibBed bed(0.042, 0.9, cfg.shift, cfg.level);
        auto& p = *bed.parametrization;
        auto model = QuantLib::ext::make_shared<ForwardMarketModel>(bed.parametrization);
        std::vector<FmmSwaptionVolTarget> targets;
        for (const auto& c : cot) {
            FmmSwaptionVolTarget t;
            t.swap = bed.coterminal(c.first);
            t.normalVol = c.second;
            t.strike = fmmForwardSwapRate(p, t.swap) - 0.003;
            std::ostringstream lbl;
            lbl << "cot_" << c.first / 4 << "y_to_5y_K-30bp";
            t.label = lbl.str();
            targets.push_back(t);
        }
        FmmSeparableVols v;
        v.segmentTimes = bed.segTimes;
        v.a = {1.0, 1.0, 1.0, 1.0};
        v.levels.assign(20, cfg.level);
        fmmSwaptionTimeDependenceBootstrap(p, v, targets, FmmSwaptionApproxMethod::EffectiveShift);
        Real worst = 0.0;
        std::ostringstream shifts;
        for (const auto& t : targets) {
            const auto r = fmmSwaptionApprox(p, t.swap, t.strike, Option::Call, FmmSwaptionApproxMethod::EffectiveShift);
            worst = std::max(worst, std::fabs(r.normalVol - t.normalVol) * 1e4);
            shifts << " " << r.effectiveShift;
        }
        BOOST_TEST_MESSAGE(cfg.label << ": analytic deal-strike fit residual " << worst << " bp; a(t) = {" << v.a[0] << ", "
                                     << v.a[1] << ", " << v.a[2] << ", " << v.a[3] << "}; effective shifts:" << shifts.str());
        BOOST_CHECK_MESSAGE(worst < 0.01, cfg.label << ": deal-strike analytic fit residual " << worst << " bp");
        BOOST_CHECK_SMALL(v.a[0] - 1.0, 1e-14);
        if (!cfg.mcLoop)
            continue;
        const auto rep = fmmMcCorrectedSwaptionBootstrap(p, model, v, targets, FmmSwaptionApproxMethod::EffectiveShift,
                                                         16384, 8, 9100, 3, 0.2);
        std::ostringstream per;
        for (const Real w : rep.worstResidualBp)
            per << " " << w;
        BOOST_TEST_MESSAGE(cfg.label << " MC-corrected loop: iterations " << rep.iterations << ", worst MC residual per iteration (bp):"
                                     << per.str() << ", converged " << rep.converged << ", runtime " << rep.runtimeSeconds
                                     << " s; a(t) = {" << v.a[0] << ", " << v.a[1] << ", " << v.a[2] << ", " << v.a[3] << "}");
        for (const auto& r : rep.rows) {
            BOOST_TEST_MESSAGE("MC-CORRECTED ROW | " << r.instrument << " | K " << r.strike << " | market " << r.marketVol * 1e4
                                                     << " bp | approx " << r.approxVol * 1e4 << " bp | mc " << r.mcVol * 1e4
                                                     << " +/- " << r.mcVolSe * 1e4 << " bp | residual " << r.residualBp
                                                     << " bp");
            BOOST_CHECK_MESSAGE(r.mcVolSe * 1e4 < 0.12, r.instrument << ": MC precision " << r.mcVolSe * 1e4 << " bp");
            BOOST_CHECK_MESSAGE(std::fabs(r.residualBp) < 3.0 * r.mcVolSe * 1e4,
                                r.instrument << ": MC residual " << r.residualBp << " bp beyond 3 s.e.");
        }
    }
}

BOOST_AUTO_TEST_CASE(testLsmEuropeanLimit) {
    BOOST_TEST_MESSAGE("A4 acceptance 1: one-right LSM 'Bermudan' matches the A3 European closed "
                       "form / approximation...");
    FmmTestBed bed(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, 3);
    FmmSwapSpec spec;
    spec.a = 8;
    spec.b = 20;
    for (Size c = 12; c <= 20; c += 4) {
        spec.fixedPayIndices.push_back(c);
        spec.fixedAccruals.push_back(1.0);
    }
    const Real K = fmmForwardSwapRate(*bed.parametrization, spec); // atm
    const auto apx = fmmSwaptionApprox(*bed.parametrization, spec, K);

    FmmCallableInstrument inst;
    inst.style = FmmCallableInstrument::Style::Enter;
    fillPayerSwapFlows(inst, *bed.parametrization, 8, 20, K);
    inst.rights.push_back({8, 8, 0.0});

    FmmLsmConfig cfg;
    cfg.valuationPaths = 262144; // sampling precision materially below 1% of the option value
    FmmLsmPricer pricer(bed.model, inst, cfg);
    const auto res = pricer.calculate();
    BOOST_TEST_MESSAGE("european limit: lsm " << res.lowerBound << " +/- " << res.lowerBoundSe << " ("
                                              << 100.0 * res.lowerBoundSe / apx.price << "% of value), approx "
                                              << apx.price << ", diff " << res.lowerBound - apx.price << " (95% CI +/- "
                                              << 1.96 * res.lowerBoundSe << "), training " << res.trainingValue
                                              << ", runtime " << res.runtimeSeconds << " s");
    BOOST_CHECK_MESSAGE(res.lowerBoundSe < 0.005 * apx.price, "sampling precision not below 0.5% of value");
    BOOST_CHECK_MESSAGE(std::fabs(res.lowerBound - apx.price) < 3.0 * res.lowerBoundSe,
                        "european limit mismatch: " << res.lowerBound - apx.price);
}

BOOST_AUTO_TEST_CASE(testLsmBermudanVsLgmGrid) {
    BOOST_TEST_MESSAGE("A4 acceptance 2: replication-mode Bermudan swaption, LSM lower bound vs "
                       "ORE's NumericLgmSwaptionEngine grid price...");

    const Date asof(19, September, 2026);
    Settings::instance().evaluationDate() = asof;
    Handle<YieldTermStructure> curve(
        QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.03, Actual365Fixed()));
    const Actual365Fixed dc;
    const Date end = asof + Period(5, Years);
    const Schedule quarterly(asof, end, Period(3, Months), NullCalendar(), Unadjusted, Unadjusted,
                             DateGeneration::Forward, false);
    const Size M = quarterly.size() - 1;
    Array rateTimes(M + 1);
    for (Size k = 0; k <= M; ++k)
        rateTimes[k] = dc.yearFraction(asof, quarterly[k]);
    Array alphaTimes(3);
    alphaTimes[0] = rateTimes[4];
    alphaTimes[1] = rateTimes[8];
    alphaTimes[2] = rateTimes[12];
    Array alpha(4);
    alpha[0] = 0.0090;
    alpha[1] = 0.0110;
    alpha[2] = 0.0100;
    alpha[3] = 0.0095;
    auto lgmParam = QuantLib::ext::make_shared<IrLgm1fPiecewiseConstantParametrization>(
        EURCurrency(), curve, alphaTimes, alpha, Array(), Array(1, 0.01));
    auto lgmModel = QuantLib::ext::make_shared<LinearGaussMarkovModel>(lgmParam);
    auto fmmParam = QuantLib::ext::make_shared<FmmParametrization>(EURCurrency(), curve, rateTimes, lgmParam);
    auto fmmModel = QuantLib::ext::make_shared<ForwardMarketModel>(fmmParam);

    // atm strike of the full 5y swap
    auto index = QuantLib::ext::make_shared<IborIndex>("FMMTEST", Period(3, Months), 0, EURCurrency(),
                                                       NullCalendar(), Unadjusted, false, dc, curve);
    const Schedule fixedSched(asof, end, Period(1, Years), NullCalendar(), Unadjusted, Unadjusted,
                              DateGeneration::Forward, false);
    const Schedule floatSched(asof, end, Period(3, Months), NullCalendar(), Unadjusted, Unadjusted,
                              DateGeneration::Forward, false);
    VanillaSwap probe(VanillaSwap::Payer, 1.0, fixedSched, 0.03, dc, floatSched, index, 0.0, dc);
    probe.setPricingEngine(QuantLib::ext::make_shared<DiscountingSwapEngine>(curve));
    const Real K = probe.fairRate();

    // LGM grid price of the Bermudan (annual exercises 1y..4y into the remaining swap)
    std::vector<Date> exDates = {quarterly[4], quarterly[8], quarterly[12], quarterly[16]};
    auto underlying = QuantLib::ext::make_shared<VanillaSwap>(VanillaSwap::Payer, 1.0, fixedSched, K, dc, floatSched,
                                                              index, 0.0, dc);
    auto swaption = QuantLib::ext::make_shared<Swaption>(
        underlying, QuantLib::ext::make_shared<BermudanExercise>(exDates));
    swaption->setPricingEngine(QuantLib::ext::make_shared<NumericLgmSwaptionEngine>(
        Handle<LinearGaussMarkovModel>(lgmModel), 7.0, 100, 7.0, 100, curve));
    const Real grid100 = swaption->NPV();
    // reference-grid convergence, reported separately from the sampling precision
    swaption->setPricingEngine(QuantLib::ext::make_shared<NumericLgmSwaptionEngine>(
        Handle<LinearGaussMarkovModel>(lgmModel), 7.0, 200, 7.0, 200, curve));
    const Real grid200 = swaption->NPV();
    swaption->setPricingEngine(QuantLib::ext::make_shared<NumericLgmSwaptionEngine>(
        Handle<LinearGaussMarkovModel>(lgmModel), 7.0, 400, 7.0, 400, curve));
    const Real grid400 = swaption->NPV();
    BOOST_TEST_MESSAGE("lgm grid convergence: 100 pts " << grid100 << ", 200 pts " << grid200 << ", 400 pts "
                                                        << grid400 << " (100 - 400: " << grid100 - grid400
                                                        << ", 200 - 400: " << grid200 - grid400 << ")");
    const Real gridPrice = grid400;

    // FMM LSM on the identical structure
    FmmCallableInstrument inst;
    inst.style = FmmCallableInstrument::Style::Enter;
    fillPayerSwapFlows(inst, *fmmParam, 0, M, K);
    // exercise into the remainder: switched flows j > settle
    inst.rights.push_back({4, 4, 0.0});
    inst.rights.push_back({8, 8, 0.0});
    inst.rights.push_back({12, 12, 0.0});
    inst.rights.push_back({16, 16, 0.0});
    // flows before the first exercise never belong to the option
    for (Size j = 1; j <= 4; ++j) {
        inst.floatWeights[j] = 0.0;
        inst.fixedFlows[j] = 0.0;
    }
    FmmLsmConfig cfg;
    cfg.trainingPaths = 32768;
    cfg.valuationPaths = 262144; // sampling precision materially below 1% of the option value
    FmmLsmPricer pricer(fmmModel, inst, cfg);
    const auto res = pricer.calculate();

    BOOST_TEST_MESSAGE("bermudan: lgm grid (400 pts) " << gridPrice << ", fmm lsm lower " << res.lowerBound << " +/- "
                                                       << res.lowerBoundSe << " (" << 100.0 * res.lowerBoundSe / gridPrice
                                                       << "% of value; diff " << res.lowerBound - gridPrice
                                                       << ", 95% CI +/- " << 1.96 * res.lowerBoundSe << "; training "
                                                       << res.trainingValue << "), runtime " << res.runtimeSeconds
                                                       << " s");
    std::ostringstream probs;
    for (Size r = 0; r < res.exerciseProbability.size(); ++r)
        probs << " " << res.exerciseProbability[r] << "+/-" << res.exerciseProbabilitySe[r];
    BOOST_TEST_MESSAGE("exercise probabilities (1y..4y):" << probs.str() << "; never "
                                                          << res.noExerciseProbability
                                                          << "; E[notice time | exercise] "
                                                          << res.expectedExerciseTime);
    // precision requirement, lower-bound property and statistical agreement (no tolerance floor)
    BOOST_CHECK_MESSAGE(res.lowerBoundSe < 0.005 * gridPrice, "sampling precision not below 0.5% of value");
    BOOST_CHECK_MESSAGE(res.lowerBound < gridPrice + 3.0 * res.lowerBoundSe,
                        "lower bound above grid price: " << res.lowerBound - gridPrice);
    BOOST_CHECK_MESSAGE(std::fabs(res.lowerBound - gridPrice) < 3.0 * res.lowerBoundSe,
                        "lsm vs grid: " << res.lowerBound - gridPrice << " beyond 3 s.e. " << 3.0 * res.lowerBoundSe);
}

BOOST_AUTO_TEST_CASE(testLsmDualBoundBermudan) {
    BOOST_TEST_MESSAGE("A4 acceptance 4: Andersen-Broadie dual upper bound and duality gap for the "
                       "replication-mode Bermudan vs the LGM grid price...");
    // same replication setup as testLsmBermudanVsLgmGrid
    const Date asof(19, September, 2026);
    Settings::instance().evaluationDate() = asof;
    Handle<YieldTermStructure> curve(
        QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.03, Actual365Fixed()));
    const Actual365Fixed dc;
    const Date end = asof + Period(5, Years);
    const Schedule quarterly(asof, end, Period(3, Months), NullCalendar(), Unadjusted, Unadjusted,
                             DateGeneration::Forward, false);
    const Size M = quarterly.size() - 1;
    Array rateTimes(M + 1);
    for (Size k = 0; k <= M; ++k)
        rateTimes[k] = dc.yearFraction(asof, quarterly[k]);
    Array alphaTimes(3);
    alphaTimes[0] = rateTimes[4];
    alphaTimes[1] = rateTimes[8];
    alphaTimes[2] = rateTimes[12];
    Array alpha(4);
    alpha[0] = 0.0090;
    alpha[1] = 0.0110;
    alpha[2] = 0.0100;
    alpha[3] = 0.0095;
    auto lgmParam = QuantLib::ext::make_shared<IrLgm1fPiecewiseConstantParametrization>(
        EURCurrency(), curve, alphaTimes, alpha, Array(), Array(1, 0.01));
    auto lgmModel = QuantLib::ext::make_shared<LinearGaussMarkovModel>(lgmParam);
    auto fmmParam = QuantLib::ext::make_shared<FmmParametrization>(EURCurrency(), curve, rateTimes, lgmParam);
    auto fmmModel = QuantLib::ext::make_shared<ForwardMarketModel>(fmmParam);
    auto index = QuantLib::ext::make_shared<IborIndex>("FMMTEST", Period(3, Months), 0, EURCurrency(),
                                                       NullCalendar(), Unadjusted, false, dc, curve);
    const Schedule fixedSched(asof, end, Period(1, Years), NullCalendar(), Unadjusted, Unadjusted,
                              DateGeneration::Forward, false);
    const Schedule floatSched(asof, end, Period(3, Months), NullCalendar(), Unadjusted, Unadjusted,
                              DateGeneration::Forward, false);
    VanillaSwap probe(VanillaSwap::Payer, 1.0, fixedSched, 0.03, dc, floatSched, index, 0.0, dc);
    probe.setPricingEngine(QuantLib::ext::make_shared<DiscountingSwapEngine>(curve));
    const Real K = probe.fairRate();
    std::vector<Date> exDates = {quarterly[4], quarterly[8], quarterly[12], quarterly[16]};
    auto underlying = QuantLib::ext::make_shared<VanillaSwap>(VanillaSwap::Payer, 1.0, fixedSched, K, dc, floatSched,
                                                              index, 0.0, dc);
    auto swaption = QuantLib::ext::make_shared<Swaption>(
        underlying, QuantLib::ext::make_shared<BermudanExercise>(exDates));
    swaption->setPricingEngine(QuantLib::ext::make_shared<NumericLgmSwaptionEngine>(
        Handle<LinearGaussMarkovModel>(lgmModel), 7.0, 100, 7.0, 100, curve));
    const Real gridPrice = swaption->NPV();

    FmmCallableInstrument inst;
    inst.style = FmmCallableInstrument::Style::Enter;
    fillPayerSwapFlows(inst, *fmmParam, 0, M, K);
    for (const Size a : {4, 8, 12, 16})
        inst.rights.push_back({a, a, 0.0});
    for (Size j = 1; j <= 4; ++j) {
        inst.floatWeights[j] = 0.0;
        inst.fixedFlows[j] = 0.0;
    }
    FmmLsmConfig cfg;
    cfg.trainingPaths = 32768;
    cfg.valuationPaths = 32768;
    FmmLsmPricer pricer(fmmModel, inst, cfg);
    const auto lsm = pricer.calculate();
    const auto dual = pricer.dualBound(1024, 128, 20260919);
    // tight upper bound = precise LSM lower bound (32k paths) + duality gap (regret-only noise)
    const Real upper = lsm.lowerBound + dual.gap;
    const Real upperSe = std::sqrt(lsm.lowerBoundSe * lsm.lowerBoundSe + dual.gapSe * dual.gapSe);
    BOOST_TEST_MESSAGE("DUALITY-GAP ROW | bermudan payer 5y annual-exercise (replication) | grid "
                       << gridPrice << " | lower(32k) " << lsm.lowerBound << " +/- " << lsm.lowerBoundSe << " | gap "
                       << dual.gap << " +/- " << dual.gapSe << " (" << 100.0 * dual.gap / gridPrice
                       << "% of value) | upper = lower + gap " << upper << " +/- " << upperSe << " | outer-sample "
                       << dual.lowerBound << " +/- " << dual.lowerBoundSe << " | outer " << dual.outerPaths
                       << " x inner " << dual.innerPaths << " | " << dual.runtimeSeconds << " s");
    // the true value (grid) must lie between the bounds within noise
    BOOST_CHECK_MESSAGE(upper > gridPrice - 3.0 * upperSe, "upper bound below grid price: " << upper - gridPrice);
    BOOST_CHECK_MESSAGE(lsm.lowerBound < gridPrice + 3.0 * lsm.lowerBoundSe,
                        "lower bound above grid price: " << lsm.lowerBound - gridPrice);
    BOOST_CHECK_MESSAGE(dual.gap > -3.0 * dual.gapSe, "negative duality gap beyond noise: " << dual.gap);
    // a policy this good should leave a small gap (reported; sanity bound 2% of value)
    BOOST_CHECK_MESSAGE(dual.gap < 0.02 * gridPrice + 3.0 * dual.gapSe,
                        "duality gap unexpectedly large: " << dual.gap);
}

BOOST_AUTO_TEST_CASE(testLsmCancellableParity) {
    BOOST_TEST_MESSAGE("A4 acceptance 3: cancellable swap priced directly equals swap + Bermudan "
                       "swaption within 3 s.e....");
    FmmTestBed bed(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, 3);
    const Real K = 0.031;
    const std::vector<FmmCallableInstrument::Right> rights = {{4, 4, 0.0}, {8, 8, 0.0}, {12, 12, 0.0}, {16, 16, 0.0}};

    // direct: payer swap (receive float, pay fixed) with cancellation rights, holder = float receiver
    FmmCallableInstrument cancellable;
    cancellable.style = FmmCallableInstrument::Style::Cancel;
    fillPayerSwapFlows(cancellable, *bed.parametrization, 0, 20, K);
    cancellable.rights = rights;
    FmmLsmPricer directPricer(bed.model, cancellable);
    const auto direct = directPricer.calculate();

    // decomposition: full payer swap + Bermudan option to ENTER the offsetting receiver remainder
    FmmCallableInstrument receiverOpt;
    receiverOpt.style = FmmCallableInstrument::Style::Enter;
    fillPayerSwapFlows(receiverOpt, *bed.parametrization, 0, 20, K, -1.0); // receiver remainder
    receiverOpt.rights = rights;
    FmmLsmConfig cfgB;
    cfgB.trainingSeed = 77;
    cfgB.valuationSeed = 787878; // independent runs
    FmmLsmPricer optPricer(bed.model, receiverOpt, cfgB);
    const auto opt = optPricer.calculate();

    const Real swapPv = direct.underlyingValue; // curve value of the full payer swap
    const Real lhs = direct.lowerBoundCv;       // control variate on the known swap value
    const Real rhs = swapPv + opt.lowerBound;
    const Real se = std::sqrt(direct.lowerBoundCvSe * direct.lowerBoundCvSe + opt.lowerBoundSe * opt.lowerBoundSe);
    BOOST_TEST_MESSAGE("cancellable direct raw " << direct.lowerBound << " +/- " << direct.lowerBoundSe << " (underlying MC - curve "
                                                 << direct.underlyingValueMc - swapPv << "), control-variate " << lhs
                                                 << " +/- " << direct.lowerBoundCvSe << "; swap " << swapPv
                                                 << " + bermudan " << opt.lowerBound << " = " << rhs << "; diff "
                                                 << lhs - rhs << " vs 3 s.e. " << 3.0 * se);
    BOOST_CHECK_MESSAGE(std::fabs(lhs - rhs) < 3.0 * se, "cancellable parity violated: " << lhs - rhs);
    pairedIdentityCheck(directPricer, optPricer, swapPv, "cancellable swap parity", opt.lowerBound);
}

namespace {

// prints one duality-gap table row (control-variate values, which equal the raw ones for Enter
// style) and checks the bound ordering
void gapRow(FmmLsmPricer& pricer, const std::string& label, const Real lsmLower, const Real lsmLowerSe) {
    const auto dual = pricer.dualBound(512, 64, 20260920);
    // tight upper bound = precise LSM lower bound + duality gap (regret-only noise)
    const Real upper = lsmLower + dual.gap;
    const Real upperSe = std::sqrt(lsmLowerSe * lsmLowerSe + dual.gapSe * dual.gapSe);
    BOOST_TEST_MESSAGE("DUALITY-GAP ROW | " << label << " | lower(32k, cv) " << lsmLower << " +/- " << lsmLowerSe
                                            << " | gap " << dual.gap << " +/- " << dual.gapSe << " ("
                                            << 100.0 * dual.gap / std::fabs(lsmLower)
                                            << "% of |value|) | upper = lower + gap " << upper << " +/- " << upperSe
                                            << " | outer-sample (cv) " << dual.lowerBoundCv << " +/- "
                                            << dual.lowerBoundCvSe << " | " << dual.runtimeSeconds << " s");
    BOOST_CHECK_MESSAGE(dual.gap > -3.0 * dual.gapSe, label << ": negative duality gap beyond noise, " << dual.gap);
}

} // namespace

BOOST_AUTO_TEST_CASE(testLsmCallableBondIdentity) {
    BOOST_TEST_MESSAGE("A4 product 3: callable fixed-rate bond (issuer, par calls) == bond - Bermudan "
                       "receiver swaption at zero issuer spread; spread sensitivity; gap row...");
    FmmTestBed bed(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, 3);
    const auto& p = *bed.parametrization;
    const Real coupon = 0.031;
    const std::vector<FmmCallableInstrument::Right> rights = {{4, 4, -1.0}, {8, 8, -1.0}, {12, 12, -1.0}, {16, 16, -1.0}};

    // issuer's callable bond: pays annual coupons and principal, may call at par (fee -1 at settle)
    FmmCallableInstrument bond;
    bond.style = FmmCallableInstrument::Style::Cancel;
    bond.lastFlowIdx = 20;
    bond.fixedFlows.assign(21, 0.0);
    bond.floatWeights.assign(21, 0.0);
    for (Size c = 4; c <= 20; c += 4)
        bond.fixedFlows[c] -= coupon * (p.rateTime(c) - p.rateTime(c - 4));
    bond.fixedFlows[20] -= 1.0;
    bond.rights = rights;
    FmmLsmConfig cfg;
    cfg.trainingPaths = 32768;
    cfg.valuationPaths = 32768;
    FmmLsmPricer bondPricer(bed.model, bond, cfg);
    const auto bondRes = bondPricer.calculate();

    // decomposition: straight bond liability + Bermudan option to ENTER the receiver remainder
    // (receive fixed coupon, pay float) - exact at zero spread since the float leg + par is par
    FmmCallableInstrument recv;
    recv.style = FmmCallableInstrument::Style::Enter;
    fillPayerSwapFlows(recv, p, 0, 20, coupon, -1.0);
    for (auto& r : recv.rights = rights)
        r.feeFlow = 0.0;
    FmmLsmConfig cfgB = cfg;
    cfgB.trainingSeed = 111;
    cfgB.valuationSeed = 222222;
    FmmLsmPricer optPricer(bed.model, recv, cfgB);
    const auto optRes = optPricer.calculate();

    const Real lhs = bondRes.lowerBoundCv; // control variate on the known straight-bond value
    const Real rhs = bondRes.underlyingValue + optRes.lowerBound;
    const Real se =
        std::sqrt(bondRes.lowerBoundCvSe * bondRes.lowerBoundCvSe + optRes.lowerBoundSe * optRes.lowerBoundSe);
    BOOST_TEST_MESSAGE("callable bond (issuer) raw " << bondRes.lowerBound << " +/- " << bondRes.lowerBoundSe
                                                      << " (underlying MC - curve "
                                                      << bondRes.underlyingValueMc - bondRes.underlyingValue
                                                      << "), control-variate " << lhs << " +/- " << bondRes.lowerBoundCvSe
                                                      << "; straight bond " << bondRes.underlyingValue
                                                      << " + bermudan receiver " << optRes.lowerBound << " = " << rhs
                                                      << "; diff " << lhs - rhs << " vs 3 s.e. " << 3.0 * se);
    BOOST_CHECK_MESSAGE(std::fabs(lhs - rhs) < 3.0 * se, "callable bond identity violated: " << lhs - rhs);
    pairedIdentityCheck(bondPricer, optPricer, bondRes.underlyingValue, "callable fixed bond", optRes.lowerBound);
    gapRow(bondPricer, "callable fixed bond 5y annual par calls (issuer)", bondRes.lowerBoundCv, bondRes.lowerBoundCvSe);

    // issuer spread 50 bp: liability PV magnitude shrinks (value less negative), option changes
    FmmCallableInstrument bondSpread = bond;
    bondSpread.issuerSpread = 0.005;
    FmmLsmPricer spreadPricer(bed.model, bondSpread, cfg);
    const auto spreadRes = spreadPricer.calculate();
    BOOST_TEST_MESSAGE("with 50 bp issuer spread: callable bond " << spreadRes.lowerBound << " (straight "
                                                                  << spreadRes.underlyingValue << ")");
    BOOST_CHECK_MESSAGE(spreadRes.underlyingValue > bondRes.underlyingValue,
                        "issuer spread must reduce the liability magnitude");
    BOOST_CHECK_MESSAGE(spreadRes.lowerBound > bondRes.lowerBound - 3.0 * (spreadRes.lowerBoundSe + bondRes.lowerBoundSe),
                        "callable value with spread not consistent");
}

BOOST_AUTO_TEST_CASE(testLsmAccretingNoteAndSwap) {
    BOOST_TEST_MESSAGE("A4 product 4: callable accreting zero note (calls at accreted value) == "
                       "zero liability - Bermudan zero-bond option; accreting cancellable swap parity; "
                       "gap rows...");
    FmmTestBed bed(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, 3);
    const auto& p = *bed.parametrization;
    const Real accRate = 0.032;
    auto accreted = [&](const Time t) { return std::pow(1.0 + accRate, t); };

    // non-call 2y, then annual calls at accreted value; single payment of accreted principal at 5y
    FmmCallableInstrument note;
    note.style = FmmCallableInstrument::Style::Cancel;
    note.lastFlowIdx = 20;
    note.fixedFlows.assign(21, 0.0);
    note.floatWeights.assign(21, 0.0);
    note.fixedFlows[20] = -accreted(p.rateTime(20));
    for (const Size a : {8, 12, 16})
        note.rights.push_back({a, a, -accreted(p.rateTime(a))});
    FmmLsmConfig cfg;
    cfg.trainingPaths = 32768;
    cfg.valuationPaths = 32768;
    FmmLsmPricer notePricer(bed.model, note, cfg);
    const auto noteRes = notePricer.calculate();

    // the call option as an Enter instrument: on exercise the issuer stops the maturity payment
    // (receives it back) and pays the accreted value at settle
    FmmCallableInstrument opt;
    opt.style = FmmCallableInstrument::Style::Enter;
    opt.lastFlowIdx = 20;
    opt.fixedFlows.assign(21, 0.0);
    opt.floatWeights.assign(21, 0.0);
    opt.fixedFlows[20] = +accreted(p.rateTime(20));
    for (const Size a : {8, 12, 16})
        opt.rights.push_back({a, a, -accreted(p.rateTime(a))});
    FmmLsmConfig cfgB = cfg;
    cfgB.trainingSeed = 333;
    cfgB.valuationSeed = 444444;
    FmmLsmPricer optPricer(bed.model, opt, cfgB);
    const auto optRes = optPricer.calculate();

    const Real lhs = noteRes.lowerBoundCv; // control variate on the known zero-liability value
    const Real rhs = noteRes.underlyingValue + optRes.lowerBound;
    const Real se =
        std::sqrt(noteRes.lowerBoundCvSe * noteRes.lowerBoundCvSe + optRes.lowerBoundSe * optRes.lowerBoundSe);
    BOOST_TEST_MESSAGE("accreting callable note raw " << noteRes.lowerBound << " +/- " << noteRes.lowerBoundSe
                                                      << ", control-variate " << lhs << " +/- " << noteRes.lowerBoundCvSe
                                                      << "; zero liability " << noteRes.underlyingValue
                                                      << " + bermudan zero-bond option " << optRes.lowerBound << " = "
                                                      << rhs << "; diff " << lhs - rhs << " vs 3 s.e. " << 3.0 * se);
    BOOST_CHECK_MESSAGE(std::fabs(lhs - rhs) < 3.0 * se, "accreting note identity violated: " << lhs - rhs);
    pairedIdentityCheck(notePricer, optPricer, noteRes.underlyingValue, "callable accreting zero note 5yNC2",
                        optRes.lowerBound);
    std::ostringstream probs;
    for (Size r = 0; r < noteRes.exerciseProbability.size(); ++r)
        probs << " " << noteRes.exerciseProbability[r];
    BOOST_TEST_MESSAGE("accreting note call probabilities (2y,3y,4y):" << probs.str() << "; never "
                                                                       << noteRes.noExerciseProbability);
    gapRow(notePricer, "callable accreting zero note 5yNC2 (issuer)", noteRes.lowerBoundCv, noteRes.lowerBoundCvSe);

    // matching accreting cancellable swap: notional accretes, fixed coupon accRate on the
    // accreted notional, float on the accreted notional; cancel rights as the note's calls
    FmmCallableInstrument accSwap;
    accSwap.style = FmmCallableInstrument::Style::Cancel;
    accSwap.lastFlowIdx = 20;
    accSwap.fixedFlows.assign(21, 0.0);
    accSwap.floatWeights.assign(21, 0.0);
    for (Size j = 1; j <= 20; ++j)
        accSwap.floatWeights[j] = accreted(p.rateTime(j - 1)); // receive float on accreted notional
    for (Size c = 4; c <= 20; c += 4)
        accSwap.fixedFlows[c] = -accRate * (p.rateTime(c) - p.rateTime(c - 4)) * accreted(p.rateTime(c - 4));
    for (const Size a : {8, 12, 16})
        accSwap.rights.push_back({a, a, 0.0});
    FmmLsmPricer accSwapPricer(bed.model, accSwap, cfg);
    const auto accSwapRes = accSwapPricer.calculate();
    // parity vs swap + Bermudan option to enter the offsetting (receiver) accreting remainder
    FmmCallableInstrument accRecv = accSwap;
    accRecv.style = FmmCallableInstrument::Style::Enter;
    for (Size j = 1; j <= 20; ++j) {
        accRecv.floatWeights[j] *= -1.0;
        accRecv.fixedFlows[j] *= -1.0;
    }
    FmmLsmPricer accRecvPricer(bed.model, accRecv, cfgB);
    const auto accRecvRes = accRecvPricer.calculate();
    const Real diffAcc = accSwapRes.lowerBoundCv - (accSwapRes.underlyingValue + accRecvRes.lowerBound);
    const Real seAcc = std::sqrt(accSwapRes.lowerBoundCvSe * accSwapRes.lowerBoundCvSe +
                                 accRecvRes.lowerBoundSe * accRecvRes.lowerBoundSe);
    BOOST_TEST_MESSAGE("accreting cancellable swap raw " << accSwapRes.lowerBound << " +/- " << accSwapRes.lowerBoundSe
                                                         << ", control-variate " << accSwapRes.lowerBoundCv << " +/- "
                                                         << accSwapRes.lowerBoundCvSe << " vs swap+option "
                                                         << accSwapRes.underlyingValue + accRecvRes.lowerBound
                                                         << "; diff " << diffAcc << " vs 3 s.e. " << 3.0 * seAcc);
    BOOST_CHECK_MESSAGE(std::fabs(diffAcc) < 3.0 * seAcc, "accreting swap parity violated: " << diffAcc);
    pairedIdentityCheck(accSwapPricer, accRecvPricer, accSwapRes.underlyingValue, "accreting cancellable swap 5yNC2",
                        accRecvRes.lowerBound);
    gapRow(accSwapPricer, "accreting cancellable swap 5yNC2", accSwapRes.lowerBoundCv, accSwapRes.lowerBoundCvSe);
}

BOOST_AUTO_TEST_CASE(testLsmStepUpNote) {
    BOOST_TEST_MESSAGE("A4 product 5: callable step-up note (issuer, par calls) == bond - Bermudan "
                       "receiver on the step-up coupon; gap row...");
    FmmTestBed bed(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, 3);
    const auto& p = *bed.parametrization;
    const std::vector<Real> coupons = {0.025, 0.0275, 0.030, 0.0325, 0.035}; // annual step-up
    const std::vector<FmmCallableInstrument::Right> rights = {{4, 4, -1.0}, {8, 8, -1.0}, {12, 12, -1.0}, {16, 16, -1.0}};

    FmmCallableInstrument note;
    note.style = FmmCallableInstrument::Style::Cancel;
    note.lastFlowIdx = 20;
    note.fixedFlows.assign(21, 0.0);
    note.floatWeights.assign(21, 0.0);
    for (Size y = 1; y <= 5; ++y)
        note.fixedFlows[4 * y] -= coupons[y - 1] * (p.rateTime(4 * y) - p.rateTime(4 * y - 4));
    note.fixedFlows[20] -= 1.0;
    note.rights = rights;
    FmmLsmConfig cfg;
    cfg.trainingPaths = 32768;
    cfg.valuationPaths = 32768;
    FmmLsmPricer notePricer(bed.model, note, cfg);
    const auto noteRes = notePricer.calculate();

    FmmCallableInstrument recv;
    recv.style = FmmCallableInstrument::Style::Enter;
    recv.lastFlowIdx = 20;
    recv.fixedFlows.assign(21, 0.0);
    recv.floatWeights.assign(21, 0.0);
    for (Size j = 1; j <= 20; ++j)
        recv.floatWeights[j] = -1.0; // pay float
    for (Size y = 1; y <= 5; ++y)
        recv.fixedFlows[4 * y] += coupons[y - 1] * (p.rateTime(4 * y) - p.rateTime(4 * y - 4)); // receive step-up
    for (auto r : rights) {
        r.feeFlow = 0.0;
        recv.rights.push_back(r);
    }
    FmmLsmConfig cfgB = cfg;
    cfgB.trainingSeed = 555;
    cfgB.valuationSeed = 666666;
    FmmLsmPricer optPricer(bed.model, recv, cfgB);
    const auto optRes = optPricer.calculate();

    const Real lhs = noteRes.lowerBoundCv; // control variate on the known straight-note value
    const Real rhs = noteRes.underlyingValue + optRes.lowerBound;
    const Real se =
        std::sqrt(noteRes.lowerBoundCvSe * noteRes.lowerBoundCvSe + optRes.lowerBoundSe * optRes.lowerBoundSe);
    BOOST_TEST_MESSAGE("step-up callable note raw " << noteRes.lowerBound << " +/- " << noteRes.lowerBoundSe
                                                    << " (underlying MC - curve "
                                                    << noteRes.underlyingValueMc - noteRes.underlyingValue
                                                    << "), control-variate " << lhs << " +/- " << noteRes.lowerBoundCvSe
                                                    << "; straight " << noteRes.underlyingValue << " + bermudan receiver "
                                                    << optRes.lowerBound << " = " << rhs << "; diff " << lhs - rhs
                                                    << " vs 3 s.e. " << 3.0 * se);
    BOOST_CHECK_MESSAGE(std::fabs(lhs - rhs) < 3.0 * se, "step-up note identity violated: " << lhs - rhs);
    pairedIdentityCheck(notePricer, optPricer, noteRes.underlyingValue, "callable step-up note", optRes.lowerBound);
    gapRow(notePricer, "callable step-up note 5y annual par calls (issuer)", noteRes.lowerBoundCv, noteRes.lowerBoundCvSe);
}

BOOST_AUTO_TEST_CASE(testExerciseTransferConfigurable) {
    BOOST_TEST_MESSAGE("A4 product 6 (owner-specified design): configurable imported-decision "
                       "scenario - issuer-optimal note call vs the dealer's frozen swap-"
                       "cancellation decisions applied to the note; notice != settle, call fee "
                       "enforced; paired CI on common valuation paths; hedge-configuration sweep; "
                       "notice-period dual bound...");
    // VALUATION PERSPECTIVE: the ISSUER of a callable fixed-rate note (all flows signed from the
    // issuer, i.e. negative coupons/principal); a POSITIVE mismatch below is value the issuer
    // LOSES by following the imported (dealer) decisions instead of the note-optimal policy.
    FmmTestBed bed(FmmParametrization::LocalVolType::DisplacedDiffusion, 0.0025, -1.0, 3);
    const auto& p = *bed.parametrization;
    const Real coupon = 0.032;
    const std::vector<FmmCallableInstrument::Right> noteRights = {
        {8, 9, -1.0}, {12, 13, -1.0}, {16, 17, -1.0}}; // notice 1 quarter before settlement; call at par

    FmmCallableInstrument note;
    note.style = FmmCallableInstrument::Style::Cancel;
    note.lastFlowIdx = 20;
    note.fixedFlows.assign(21, 0.0);
    note.floatWeights.assign(21, 0.0);
    for (Size c = 4; c <= 20; c += 4)
        note.fixedFlows[c] -= coupon * (p.rateTime(c) - p.rateTime(c - 4));
    note.fixedFlows[20] -= 1.0; // principal
    note.rights = noteRights;

    FmmLsmPricer notePricer(bed.model, note);
    const auto noteOwn = notePricer.calculate();
    std::ostringstream ownProb;
    for (Size r = 0; r < noteOwn.exerciseProbability.size(); ++r)
        ownProb << " " << noteOwn.exerciseProbability[r];
    BOOST_TEST_MESSAGE("note own-policy " << noteOwn.lowerBound << " +/- " << noteOwn.lowerBoundSe
                                          << "; exercise probs:" << ownProb.str() << " (never "
                                          << noteOwn.noExerciseProbability << ")");

    // notice-period dual bound: rights with settle > notice use the adapted (notice-date
    // conditional expectation) exercise payoff
    {
        const auto dual = notePricer.dualBound(512, 64, 20260921);
        const Real upper = noteOwn.lowerBoundCv + dual.gap;
        const Real upperSe = std::sqrt(noteOwn.lowerBoundCvSe * noteOwn.lowerBoundCvSe + dual.gapSe * dual.gapSe);
        BOOST_TEST_MESSAGE("DUALITY-GAP ROW | callable note 5yNC2, 1q notice period (issuer) | lower(16k, cv) "
                           << noteOwn.lowerBoundCv << " +/- " << noteOwn.lowerBoundCvSe << " | gap " << dual.gap
                           << " +/- " << dual.gapSe << " | upper = lower + gap " << upper << " +/- " << upperSe
                           << " | " << dual.runtimeSeconds << " s");
        BOOST_CHECK_MESSAGE(dual.gap > -3.0 * dual.gapSe,
                            "notice-period: negative duality gap beyond noise: " << dual.gap);
        BOOST_CHECK_MESSAGE(dual.gap < 0.005 + 3.0 * dual.gapSe,
                            "notice-period duality gap unexpectedly large: " << dual.gap);
    }

    // the dealer's mirror cancellable swap in several hedge CONFIGURATIONS: aligned pairing (the
    // dealer PAYS fixed and receives float - the issuer's note-to-floating hedge - so its
    // cancellation trigger, rates falling, aligns with the note call) at the note coupon and at
    // +/- 50 bp, and the opposite pairing (the dealer receives fixed). The dealer decides on its
    // OWN swap; the issuer follows those decisions on the note. Each row is one configuration of
    // the scenario, none is a claim of universal practice.
    struct Config {
        std::string label;
        Real fixedRate;
        Real sign;
    };
    const std::vector<Config> configs = {{"aligned, dealer pays fixed = coupon", coupon, 1.0},
                                         {"aligned, dealer pays fixed = coupon - 50bp", coupon - 0.005, 1.0},
                                         {"aligned, dealer pays fixed = coupon + 50bp", coupon + 0.005, 1.0},
                                         {"opposite, dealer receives fixed = coupon", coupon, -1.0}};
    Size cfgNo = 0;
    for (const auto& c : configs) {
        FmmCallableInstrument dealerSwap;
        dealerSwap.style = FmmCallableInstrument::Style::Cancel;
        fillPayerSwapFlows(dealerSwap, p, 0, 20, c.fixedRate, c.sign);
        dealerSwap.rights = {{8, 9, 0.0}, {12, 13, 0.0}, {16, 17, 0.0}};
        FmmLsmConfig cfgD;
        cfgD.trainingSeed = static_cast<BigNatural>(91 + cfgNo++);
        FmmLsmPricer dealerPricer(bed.model, dealerSwap, cfgD);
        const auto dealerRes = dealerPricer.calculate();
        const auto noteImported = notePricer.valueWithImportedDecisions(dealerPricer, 555555);
        const auto paired = notePricer.pairedImportedDecisionDifference(dealerPricer, 666666);
        std::ostringstream dProb, impProb;
        for (Size r = 0; r < dealerRes.exerciseProbability.size(); ++r)
            dProb << " " << dealerRes.exerciseProbability[r];
        for (Size r = 0; r < noteImported.exerciseProbability.size(); ++r)
            impProb << " " << noteImported.exerciseProbability[r];
        BOOST_TEST_MESSAGE("EXERCISE-TRANSFER ROW | " << c.label << " | dealer cancel probs:" << dProb.str()
                                                      << " (never " << dealerRes.noExerciseProbability
                                                      << ") | note under imported decisions " << noteImported.lowerBound
                                                      << " +/- " << noteImported.lowerBoundSe << " (probs:" << impProb.str()
                                                      << ") | mismatch own - imported (paired) " << paired.first
                                                      << " +/- " << paired.second);
        // own-optimal must not be worse than the imported decisions beyond estimation noise
        BOOST_CHECK_MESSAGE(paired.first > -3.0 * paired.second - 2e-5,
                            c.label << ": own policy worse than imported beyond noise: " << paired.first);
        if (cfgNo == 2) {
            // the off-strike configuration must produce genuinely different exercise behaviour,
            // else the scenario is vacuous (guards against a degenerate imported policy); at the
            // note's own strike the dealer holds the same option and decides almost identically
            Real probGap = std::fabs(noteOwn.noExerciseProbability - noteImported.noExerciseProbability);
            for (Size r = 0; r < noteOwn.exerciseProbability.size(); ++r)
                probGap = std::max(probGap,
                                   std::fabs(noteOwn.exerciseProbability[r] - noteImported.exerciseProbability[r]));
            BOOST_CHECK_MESSAGE(probGap > 0.01,
                                "imported decisions indistinguishable from own policy (gap " << probGap << ")");
        }
        if (c.sign < 0.0) {
            // opposite pairing: the dealer's cancellation is in the money when rates RISE, the
            // opposite of the note call; following it must cost the issuer significantly
            BOOST_CHECK_MESSAGE(paired.first > 3.0 * paired.second,
                                "opposite configuration should show a significant positive mismatch: "
                                    << paired.first);
        }
    }
    // joint package optimization is intentionally NOT performed here (kept distinct per owner)
}

BOOST_AUTO_TEST_CASE(testLsmCallableBondVsOreLgmEngine) {
    BOOST_TEST_MESSAGE("A4 outstanding item 3: callable fixed-rate bond - FMM (replication mode) LSM "
                       "and AB dual estimates vs ORE's NumericLgmCallableBondEngine on the same LGM...");
    ReplicationBed rb;
    const auto& p = *rb.fmmParam;
    const Real coupon = 0.031;

    // ORE side (holder's perspective): annual Act/365F coupons, par redemption added by the Bond
    // base class, issuer calls at par on the 1y..4y coupon dates. In the ORE engine a coupon paid
    // on a call date is settled before the call and the accrual on a coupon date is zero, so a
    // dirty call price of 1.0 is exactly the FMM structure below (coupon at settle kept, fee -1).
    const Schedule annual(rb.asof, rb.end, Period(1, Years), NullCalendar(), Unadjusted, Unadjusted,
                          DateGeneration::Forward, false);
    const Leg coupons = FixedRateLeg(annual).withNotionals(1.0).withCouponRates(coupon, rb.dc);
    std::vector<CallableBond::CallabilityData> calls;
    for (Size y = 1; y <= 4; ++y)
        calls.push_back({annual[y], CallableBond::CallabilityData::ExerciseType::OnThisDate, 1.0,
                         CallableBond::CallabilityData::PriceType::Dirty, true});
    auto oreBond = QuantLib::ext::make_shared<CallableBond>(0, NullCalendar(), rb.asof, coupons, calls);
    // FD-solver engine: its time grid has points before the first coupon date. The convolution-
    // solver engine (timeStepsPerYear 0, event-date grid only) adds a coupon to the underlying only
    // at a grid time strictly before its accrual end, so with no grid point between today and the
    // first event it omits the coupon paying on the first event date - reproduced below for the
    // record (upstream observation, not used for the comparison).
    oreBond->setPricingEngine(QuantLib::ext::make_shared<NumericLgmCallableBondEngine>(
        Handle<LGM>(rb.lgmModel), 50.0, FdmSchemeDesc::Douglas(), 256, 48, 1e-4, 24, rb.curve));
    const Real oreNpv = oreBond->NPV();
    const Real oreStraight = oreBond->result<Real>("strippedBondNpv");
    const Real oreOption = oreBond->result<Real>("callPutValue");
    BOOST_TEST_MESSAGE("ORE NumericLgmCallableBondEngine (FD 256 x 48/y): callable bond (holder) "
                       << oreNpv << ", straight bond " << oreStraight << ", embedded call value " << oreOption);
    {
        auto conv = QuantLib::ext::make_shared<CallableBond>(0, NullCalendar(), rb.asof, coupons, calls);
        conv->setPricingEngine(QuantLib::ext::make_shared<NumericLgmCallableBondEngine>(
            Handle<LGM>(rb.lgmModel), 7.0, 100, 7.0, 100, 24, rb.curve));
        const Real convStraight = conv->result<Real>("strippedBondNpv");
        const Real firstCouponPv = coupons.front()->amount() * rb.curve->discount(coupons.front()->date());
        BOOST_TEST_MESSAGE("ORE convolution-solver engine (event-date grid): straight bond "
                           << convStraight << ", short of the FD value by " << oreStraight - convStraight
                           << " = discounted first coupon " << firstCouponPv << " (upstream observation); call value "
                           << conv->result<Real>("callPutValue"));
    }

    // FMM side (issuer's perspective, flows negative): same coupons, par call fee at the call date
    FmmCallableInstrument bond;
    bond.style = FmmCallableInstrument::Style::Cancel;
    bond.lastFlowIdx = rb.M;
    bond.fixedFlows.assign(rb.M + 1, 0.0);
    bond.floatWeights.assign(rb.M + 1, 0.0);
    for (Size c = 4; c <= rb.M; c += 4)
        bond.fixedFlows[c] -= coupon * (p.rateTime(c) - p.rateTime(c - 4));
    bond.fixedFlows[rb.M] -= 1.0;
    for (const Size a : {4, 8, 12, 16})
        bond.rights.push_back({a, a, -1.0});
    FmmLsmConfig cfg;
    cfg.trainingPaths = 32768;
    cfg.valuationPaths = 32768;
    FmmLsmPricer pricer(rb.fmmModel, bond, cfg);
    const auto lsm = pricer.calculate();
    const auto dual = pricer.dualBound(1024, 128, 20260923);
    // control-variate values (the straight bond is known exactly); embedded call value from the
    // issuer's side: LSM = lower estimate of the option, LSM + duality gap = upper estimate
    const Real fmmOptionLsm = lsm.lowerBoundCv - lsm.underlyingValue;
    const Real fmmOptionDual = fmmOptionLsm + dual.gap;
    const Real dualSe = std::sqrt(lsm.lowerBoundCvSe * lsm.lowerBoundCvSe + dual.gapSe * dual.gapSe);
    BOOST_TEST_MESSAGE("FMM replication: straight bond (curve) "
                       << -lsm.underlyingValue << "; callable bond (holder) from LSM " << -lsm.lowerBoundCv << " +/- "
                       << lsm.lowerBoundCvSe << " (raw " << -lsm.lowerBound << " +/- " << lsm.lowerBoundSe
                       << "); embedded call: LSM " << fmmOptionLsm << " +/- " << lsm.lowerBoundCvSe << ", + gap "
                       << dual.gap << " +/- " << dual.gapSe << " = " << fmmOptionDual << " vs ORE " << oreOption
                       << ": LSM diff " << fmmOptionLsm - oreOption << " ("
                       << 100.0 * (fmmOptionLsm - oreOption) / oreOption << "% of the call value), upper diff "
                       << fmmOptionDual - oreOption << "; " << lsm.runtimeSeconds + dual.runtimeSeconds << " s");
    BOOST_CHECK_MESSAGE(std::fabs(-lsm.underlyingValue - oreStraight) < 5e-5,
                        "straight bond mismatch: " << -lsm.underlyingValue - oreStraight);
    // ORE's value must lie between the FMM estimates within noise
    BOOST_CHECK_MESSAGE(oreOption > fmmOptionLsm - 3.0 * lsm.lowerBoundCvSe,
                        "ORE call value below the FMM LSM (lower) estimate: " << oreOption - fmmOptionLsm);
    BOOST_CHECK_MESSAGE(oreOption < fmmOptionDual + 3.0 * dualSe,
                        "ORE call value above the FMM dual (upper) estimate: " << oreOption - fmmOptionDual);
    // statistical agreement of the LSM estimate (no tolerance floor)
    BOOST_CHECK_MESSAGE(std::fabs(fmmOptionLsm - oreOption) < 3.0 * lsm.lowerBoundCvSe,
                        "FMM LSM vs ORE grid: " << fmmOptionLsm - oreOption << " beyond 3 s.e.");

    // issuer spread 50 bp on both sides: ORE's discounting spread (zero spread on the effective
    // discount curve) vs the FMM deterministic issuer spread
    {
        const Real s = 0.005;
        auto oreSpread = QuantLib::ext::make_shared<CallableBond>(0, NullCalendar(), rb.asof, coupons, calls);
        oreSpread->setPricingEngine(QuantLib::ext::make_shared<NumericLgmCallableBondEngine>(
            Handle<LGM>(rb.lgmModel), 50.0, FdmSchemeDesc::Douglas(), 256, 48, 1e-4, 24, rb.curve,
            Handle<Quote>(QuantLib::ext::make_shared<SimpleQuote>(s))));
        const Real oreS = oreSpread->NPV();
        const Real oreSStraight = oreSpread->result<Real>("strippedBondNpv");
        const Real oreSOption = oreSpread->result<Real>("callPutValue");
        FmmCallableInstrument bondS = bond;
        bondS.issuerSpread = s;
        FmmLsmPricer pricerS(rb.fmmModel, bondS, cfg);
        const auto lsmS = pricerS.calculate();
        const auto dualS = pricerS.dualBound(512, 64, 20260924);
        const Real optS = lsmS.lowerBoundCv - lsmS.underlyingValue;
        const Real optSUpper = optS + dualS.gap;
        const Real seS = std::sqrt(lsmS.lowerBoundCvSe * lsmS.lowerBoundCvSe + dualS.gapSe * dualS.gapSe);
        BOOST_TEST_MESSAGE("issuer spread 50 bp: ORE callable " << oreS << ", straight " << oreSStraight << ", call "
                                                                << oreSOption << " | FMM straight " << -lsmS.underlyingValue
                                                                << ", call LSM " << optS << " +/- " << lsmS.lowerBoundCvSe
                                                                << ", + gap " << dualS.gap << " = " << optSUpper << " +/- "
                                                                << seS);
        BOOST_CHECK_MESSAGE(std::fabs(-lsmS.underlyingValue - oreSStraight) < 5e-5,
                            "straight bond with spread mismatch: " << -lsmS.underlyingValue - oreSStraight);
        BOOST_CHECK_MESSAGE(oreSOption > optS - 3.0 * lsmS.lowerBoundCvSe && oreSOption < optSUpper + 3.0 * seS,
                            "ORE call value with spread outside the FMM bounds: " << oreSOption - optS);
    }
}

BOOST_AUTO_TEST_CASE(testMultiFactorVsLgmCalibratedBasket) {
    BOOST_TEST_MESSAGE("A4 acceptance 5: 3-factor FMM vs 1F LGM, both calibrated to the SAME "
                       "coterminal basket (real USD-SOFR ATM vols); differences reported, not tuned...");

    const Date asof(19, September, 2026);
    Settings::instance().evaluationDate() = asof;
    Handle<YieldTermStructure> curve(
        QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.03, Actual365Fixed()));
    const Actual365Fixed dc;
    const Date end = asof + Period(5, Years);
    const Schedule quarterly(asof, end, Period(3, Months), NullCalendar(), Unadjusted, Unadjusted,
                             DateGeneration::Forward, false);
    const Size M = quarterly.size() - 1;
    Array rateTimes(M + 1);
    for (Size k = 0; k <= M; ++k)
        rateTimes[k] = dc.yearFraction(asof, quarterly[k]);
    auto index = QuantLib::ext::make_shared<IborIndex>("FMMTEST", Period(3, Months), 0, EURCurrency(),
                                                       NullCalendar(), Unadjusted, false, dc, curve);

    // shared basket: coterminal-to-5y ATM normal vols (Products marketdata.csv, 2025-02-10)
    const std::vector<std::pair<Size, Real>> basket = {
        {4, 0.01029015}, {8, 0.0105039}, {12, 0.0105350}, {16, 0.0104501}};

    // ---- LGM: piecewise-constant alpha at the basket expiries, kappa 0.01, ORE helpers ----
    Array alphaTimes(3);
    alphaTimes[0] = rateTimes[4];
    alphaTimes[1] = rateTimes[8];
    alphaTimes[2] = rateTimes[12];
    auto lgmParam = QuantLib::ext::make_shared<IrLgm1fPiecewiseConstantParametrization>(
        EURCurrency(), curve, alphaTimes, Array(4, 0.01), Array(), Array(1, 0.01));
    auto lgmModel = QuantLib::ext::make_shared<LinearGaussMarkovModel>(lgmParam);
    auto lgmEngine = QuantLib::ext::make_shared<AnalyticLgmSwaptionEngine>(lgmParam, curve);
    std::vector<QuantLib::ext::shared_ptr<BlackCalibrationHelper>> helpers;
    for (const auto& b : basket) {
        auto h = QuantLib::ext::make_shared<SwaptionHelper>(
            Period(static_cast<Integer>(b.first / 4), Years), Period(static_cast<Integer>(5 - b.first / 4), Years),
            Handle<Quote>(QuantLib::ext::make_shared<SimpleQuote>(b.second)), index, Period(1, Years), dc, dc,
            curve, BlackCalibrationHelper::RelativePriceError, Null<Real>(), 1.0, Normal, 0.0);
        h->setPricingEngine(lgmEngine);
        helpers.push_back(h);
    }
    LevenbergMarquardt lm;
    EndCriteria ec(1000, 100, 1e-8, 1e-8, 1e-8);
    lgmModel->calibrateVolatilitiesIterative(helpers, lm, ec);
    Real worstLgmBp = 0.0;
    for (Size i = 0; i < helpers.size(); ++i) {
        const Real mv = helpers[i]->marketValue(), model = helpers[i]->modelValue();
        const Real relErr = (model - mv) / mv;
        BOOST_TEST_MESSAGE("LGM calib helper " << i << ": market " << mv << " model " << model << " rel err "
                                                << relErr);
        worstLgmBp = std::max(worstLgmBp, std::fabs(relErr));
    }
    BOOST_CHECK_MESSAGE(worstLgmBp < 1e-3, "LGM calibration did not converge: worst rel err " << worstLgmBp);

    // ---- FMM 3F: strategy (b) to the same basket (exact), correlation rhoInf 0.6 / beta 0.08 ----
    Array shifts(M);
    for (Size j = 1; j <= M; ++j)
        shifts[j - 1] = 1.0 / (rateTimes[j] - rateTimes[j - 1]);
    std::vector<Array> levels(M, Array(4, 0.0025));
    auto fmm3f = QuantLib::ext::make_shared<FmmParametrization>(EURCurrency(), curve, rateTimes, shifts, alphaTimes,
                                                                levels,
                                                                FmmParametrization::LocalVolType::DisplacedDiffusion,
                                                                0.6, 0.08, 3);
    std::vector<FmmSwaptionVolTarget> targets;
    for (const auto& b : basket) {
        FmmSwaptionVolTarget t;
        t.swap.a = b.first;
        t.swap.b = M;
        for (Size c = b.first + 4; c <= M; c += 4) {
            t.swap.fixedPayIndices.push_back(c);
            t.swap.fixedAccruals.push_back(rateTimes[c] - rateTimes[c - 4]);
        }
        t.normalVol = b.second;
        targets.push_back(t);
    }
    FmmSeparableVols v;
    v.segmentTimes = alphaTimes;
    v.a = {1.0, 1.0, 1.0, 1.0};
    v.levels.assign(M, 0.0025);
    fmmSwaptionTimeDependenceBootstrap(*fmm3f, v, targets);
    Real worstFmmBp = 0.0;
    for (const auto& t : targets)
        worstFmmBp = std::max(
            worstFmmBp,
            std::fabs(fmmSwaptionApprox(*fmm3f, t.swap, fmmForwardSwapRate(*fmm3f, t.swap)).normalVol - t.normalVol) * 1e4);
    BOOST_TEST_MESSAGE("FMM 3F calibration worst residual " << worstFmmBp << " bp; a(t) = {" << v.a[0] << ", " << v.a[1]
                                                           << ", " << v.a[2] << ", " << v.a[3] << "}");
    BOOST_CHECK_MESSAGE(worstFmmBp < 0.01, "FMM calibration residual " << worstFmmBp << " bp");

    // ---- the benchmark trade: Bermudan payer, ATM on the 5y swap, annual exercise ----
    const Schedule fixedSched(asof, end, Period(1, Years), NullCalendar(), Unadjusted, Unadjusted,
                              DateGeneration::Forward, false);
    const Schedule floatSched(asof, end, Period(3, Months), NullCalendar(), Unadjusted, Unadjusted,
                              DateGeneration::Forward, false);
    VanillaSwap probe(VanillaSwap::Payer, 1.0, fixedSched, 0.03, dc, floatSched, index, 0.0, dc);
    probe.setPricingEngine(QuantLib::ext::make_shared<DiscountingSwapEngine>(curve));
    const Real K = probe.fairRate();
    std::vector<Date> exDates = {quarterly[4], quarterly[8], quarterly[12], quarterly[16]};
    auto underlying = QuantLib::ext::make_shared<VanillaSwap>(VanillaSwap::Payer, 1.0, fixedSched, K, dc, floatSched,
                                                              index, 0.0, dc);
    auto swaption = QuantLib::ext::make_shared<Swaption>(
        underlying, QuantLib::ext::make_shared<BermudanExercise>(exDates));
    swaption->setPricingEngine(QuantLib::ext::make_shared<NumericLgmSwaptionEngine>(
        Handle<LinearGaussMarkovModel>(lgmModel), 7.0, 100, 7.0, 100, curve));
    const Real lgmBermudan = swaption->NPV();

    auto makeInst = [&](const FmmParametrization& p) {
        FmmCallableInstrument inst;
        inst.style = FmmCallableInstrument::Style::Enter;
        fillPayerSwapFlows(inst, p, 0, M, K);
        for (const Size a : {4, 8, 12, 16})
            inst.rights.push_back({a, a, 0.0});
        for (Size j = 1; j <= 4; ++j) {
            inst.floatWeights[j] = 0.0;
            inst.fixedFlows[j] = 0.0;
        }
        return inst;
    };
    FmmLsmConfig cfg;
    cfg.trainingPaths = 32768;
    cfg.valuationPaths = 32768;

    // replication-mode FMM of the CALIBRATED LGM: the attribution anchor (must match the grid)
    auto fmmRep = QuantLib::ext::make_shared<FmmParametrization>(EURCurrency(), curve, rateTimes, lgmParam);
    auto repModel = QuantLib::ext::make_shared<ForwardMarketModel>(fmmRep);
    FmmLsmPricer repPricer(repModel, makeInst(*fmmRep), cfg);
    const auto rep = repPricer.calculate();

    // multi-factor FMM on the same basket
    auto model3f = QuantLib::ext::make_shared<ForwardMarketModel>(fmm3f);
    FmmLsmPricer pricer3f(model3f, makeInst(*fmm3f), cfg);
    const auto mf = pricer3f.calculate();
    const auto dual3f = pricer3f.dualBound(512, 64, 20260921);
    // paired common random numbers: both models driven by the same Gaussian draws (seed 424242);
    // the pathwise 3F - replication difference has far less noise than the unpaired difference
    std::pair<Real, Real> paired3f;
    {
        std::vector<Real> v3, t3, vr, tr;
        pricer3f.pathValues(cfg.valuationSeed, v3, t3);
        repPricer.pathValues(cfg.valuationSeed, vr, tr);
        IncrementalStatistics d;
        for (Size n = 0; n < v3.size(); ++n)
            d.add(v3[n] - vr[n]);
        paired3f = {d.mean(), d.errorEstimate()};
    }
    BOOST_TEST_MESSAGE("PAIRED-CRN CHALLENGER | 3F - replication on common draws: "
                       << paired3f.first << " +/- " << paired3f.second << " (" << 100.0 * paired3f.first / lgmBermudan
                       << "% of value, 95% CI +/- " << 100.0 * 1.96 * paired3f.second / lgmBermudan << "%)");

    BOOST_TEST_MESSAGE("CHALLENGER ROW | bermudan payer 5y annual, shared coterminal basket | LGM grid "
                       << lgmBermudan << " | FMM-replication LSM " << rep.lowerBound << " +/- " << rep.lowerBoundSe
                       << " | FMM-3F LSM lower " << mf.lowerBound << " +/- " << mf.lowerBoundSe << ", duality gap "
                       << dual3f.gap << " +/- " << dual3f.gapSe << ", upper = lower + gap "
                       << mf.lowerBound + dual3f.gap << " +/- "
                       << std::sqrt(mf.lowerBoundSe * mf.lowerBoundSe + dual3f.gapSe * dual3f.gapSe));
    BOOST_TEST_MESSAGE("3F minus LGM: " << mf.lowerBound - lgmBermudan << " (" << 100.0 * (mf.lowerBound - lgmBermudan) / lgmBermudan
                                        << "% of value); exercise probs 3F:" << mf.exerciseProbability[0] << " "
                                        << mf.exerciseProbability[1] << " " << mf.exerciseProbability[2] << " "
                                        << mf.exerciseProbability[3] << " vs replication:" << rep.exerciseProbability[0]
                                        << " " << rep.exerciseProbability[1] << " " << rep.exerciseProbability[2] << " "
                                        << rep.exerciseProbability[3]);
    BOOST_TEST_MESSAGE("attribution: both models reprice the same 4 ATM coterminals (LGM rel err < 1e-3, FMM < 0.01 bp); "
                       "the residual difference is the multi-factor (decorrelation) effect plus the DD-vs-Gaussian "
                       "distribution of forward rates at the shift 1/tau, not implementation noise "
                       "(replication anchor vs grid: " << rep.lowerBound - lgmBermudan << ")");
    // the anchor must reproduce the LGM grid within the LSM tolerance
    BOOST_CHECK_MESSAGE(std::fabs(rep.lowerBound - lgmBermudan) < std::max(3.0 * rep.lowerBoundSe, 0.01 * lgmBermudan),
                        "replication anchor off the LGM grid: " << rep.lowerBound - lgmBermudan);
    // no assertion on the 3F difference: it is a reported finding by design (do not tune away)
    BOOST_CHECK(dual3f.upperBound > mf.lowerBound - 3.0 * (dual3f.upperBoundSe + mf.lowerBoundSe));
}

BOOST_AUTO_TEST_CASE(testLsmDualBoundProductionBudget) {
    // minutes-long: opt in with FMM_LONG_TESTS=1; produces the production-budget rows of the
    // duality-gap table (inner-path convergence at fixed outer budget)
    if (std::getenv("FMM_LONG_TESTS") == nullptr) {
        BOOST_TEST_MESSAGE("A4 duality-gap table, production budget: SKIPPED (set FMM_LONG_TESTS=1 to run)");
        BOOST_CHECK(true);
        return;
    }
    BOOST_TEST_MESSAGE("A4 duality-gap table, production inner-path budget...");
    ReplicationBed rb;
    const auto& p = *rb.fmmParam;
    const std::vector<std::pair<Size, Size>> budgets = {{2048, 256}, {2048, 1024}};
    auto rows = [&](FmmLsmPricer& pricer, const std::string& label, const Real lsmLower, const Real lsmLowerSe,
                    const Real truth) {
        for (const auto& b : budgets) {
            const auto dual = pricer.dualBound(b.first, b.second, 20260922);
            const Real upper = lsmLower + dual.gap;
            const Real upperSe = std::sqrt(lsmLowerSe * lsmLowerSe + dual.gapSe * dual.gapSe);
            BOOST_TEST_MESSAGE("DUALITY-GAP ROW (production) | "
                               << label << " | lower(32k, cv) " << lsmLower << " +/- " << lsmLowerSe << " | gap "
                               << dual.gap << " +/- " << dual.gapSe << " | upper = lower + gap " << upper << " +/- "
                               << upperSe << " | truth " << (truth != Null<Real>() ? std::to_string(truth) : "n/a")
                               << " | outer " << b.first << " x inner " << b.second << " | " << dual.runtimeSeconds
                               << " s");
            BOOST_CHECK_MESSAGE(dual.gap > -3.0 * dual.gapSe,
                                label << ": negative duality gap beyond noise: " << dual.gap);
            if (truth != Null<Real>())
                BOOST_CHECK_MESSAGE(upper > truth - 3.0 * upperSe,
                                    label << ": upper bound below the reference value: " << upper - truth);
        }
    };

    // (a) replication-mode Bermudan payer, LGM grid price as truth
    auto index = QuantLib::ext::make_shared<IborIndex>("FMMTEST", Period(3, Months), 0, EURCurrency(), NullCalendar(),
                                                       Unadjusted, false, rb.dc, rb.curve);
    const Schedule fixedSched(rb.asof, rb.end, Period(1, Years), NullCalendar(), Unadjusted, Unadjusted,
                              DateGeneration::Forward, false);
    const Schedule floatSched(rb.asof, rb.end, Period(3, Months), NullCalendar(), Unadjusted, Unadjusted,
                              DateGeneration::Forward, false);
    VanillaSwap probe(VanillaSwap::Payer, 1.0, fixedSched, 0.03, rb.dc, floatSched, index, 0.0, rb.dc);
    probe.setPricingEngine(QuantLib::ext::make_shared<DiscountingSwapEngine>(rb.curve));
    const Real K = probe.fairRate();
    std::vector<Date> exDates = {rb.quarterly[4], rb.quarterly[8], rb.quarterly[12], rb.quarterly[16]};
    auto underlying = QuantLib::ext::make_shared<VanillaSwap>(VanillaSwap::Payer, 1.0, fixedSched, K, rb.dc,
                                                              floatSched, index, 0.0, rb.dc);
    auto swaption =
        QuantLib::ext::make_shared<Swaption>(underlying, QuantLib::ext::make_shared<BermudanExercise>(exDates));
    swaption->setPricingEngine(QuantLib::ext::make_shared<NumericLgmSwaptionEngine>(
        Handle<LinearGaussMarkovModel>(rb.lgmModel), 7.0, 100, 7.0, 100, rb.curve));
    const Real gridPrice = swaption->NPV();
    FmmCallableInstrument inst;
    inst.style = FmmCallableInstrument::Style::Enter;
    fillPayerSwapFlows(inst, p, 0, rb.M, K);
    for (const Size a : {4, 8, 12, 16})
        inst.rights.push_back({a, a, 0.0});
    for (Size j = 1; j <= 4; ++j) {
        inst.floatWeights[j] = 0.0;
        inst.fixedFlows[j] = 0.0;
    }
    FmmLsmConfig cfg;
    cfg.trainingPaths = 32768;
    cfg.valuationPaths = 32768;
    FmmLsmPricer pricer(rb.fmmModel, inst, cfg);
    const auto lsm = pricer.calculate();
    rows(pricer, "bermudan payer 5y annual-exercise (replication)", lsm.lowerBound, lsm.lowerBoundSe, gridPrice);

    // (b) callable fixed-rate bond (issuer, annual par calls), no independent truth
    FmmCallableInstrument bond;
    bond.style = FmmCallableInstrument::Style::Cancel;
    bond.lastFlowIdx = rb.M;
    bond.fixedFlows.assign(rb.M + 1, 0.0);
    bond.floatWeights.assign(rb.M + 1, 0.0);
    for (Size c = 4; c <= rb.M; c += 4)
        bond.fixedFlows[c] -= 0.031 * (p.rateTime(c) - p.rateTime(c - 4));
    bond.fixedFlows[rb.M] -= 1.0;
    for (const Size a : {4, 8, 12, 16})
        bond.rights.push_back({a, a, -1.0});
    FmmLsmPricer bondPricer(rb.fmmModel, bond, cfg);
    const auto bondRes = bondPricer.calculate();
    rows(bondPricer, "callable fixed bond 5y annual par calls (issuer, replication)", bondRes.lowerBoundCv,
         bondRes.lowerBoundCvSe, Null<Real>());
}

BOOST_AUTO_TEST_CASE(testLsmAccretingNote30NC5) {
    // minutes-long (30y semiannual grid, 25 annual rights): opt in with FMM_LONG_TESTS=1
    if (std::getenv("FMM_LONG_TESTS") == nullptr) {
        BOOST_TEST_MESSAGE("A4 product 4, 30Y non-call-5 accreting example: SKIPPED (set FMM_LONG_TESTS=1 to run)");
        BOOST_CHECK(true);
        return;
    }
    BOOST_TEST_MESSAGE("A4 product 4 (revised plan): 30Y non-call-5 callable accreting zero note, annual calls "
                       "at accreted value; identity vs zero liability - Bermudan zero-bond option (independent "
                       "and paired); matching accreting cancellable swap; gap row...");
    Settings::instance().evaluationDate() = Date(19, September, 2026);
    Handle<YieldTermStructure> curve(
        QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.03, Actual365Fixed()));
    const Size M = 60; // semiannual grid to 30y
    Array rateTimes(M + 1);
    for (Size k = 0; k <= M; ++k)
        rateTimes[k] = 0.5 * static_cast<Real>(k);
    Array shifts(M, 2.0); // 1/tau
    std::vector<Array> volLevels(M, Array(1, 0.0025));
    auto param = QuantLib::ext::make_shared<FmmParametrization>(EURCurrency(), curve, rateTimes, shifts, Array(),
                                                                volLevels,
                                                                FmmParametrization::LocalVolType::DisplacedDiffusion,
                                                                0.6, 0.08, 3);
    auto model = QuantLib::ext::make_shared<ForwardMarketModel>(param);
    const auto& p = *param;
    const Real accRate = 0.032;
    auto accreted = [&](const Time t) { return std::pow(1.0 + accRate, t); };
    std::vector<Size> callIdx;
    for (Size a = 10; a <= 58; a += 2)
        callIdx.push_back(a); // 5y .. 29y annually

    FmmCallableInstrument note;
    note.style = FmmCallableInstrument::Style::Cancel;
    note.lastFlowIdx = M;
    note.fixedFlows.assign(M + 1, 0.0);
    note.floatWeights.assign(M + 1, 0.0);
    note.fixedFlows[M] = -accreted(p.rateTime(M));
    for (const Size a : callIdx)
        note.rights.push_back({a, a, -accreted(p.rateTime(a))});
    FmmLsmConfig cfg;
    cfg.trainingPaths = 16384;
    cfg.valuationPaths = 16384;
    FmmLsmPricer notePricer(model, note, cfg);
    const auto noteRes = notePricer.calculate();

    FmmCallableInstrument opt;
    opt.style = FmmCallableInstrument::Style::Enter;
    opt.lastFlowIdx = M;
    opt.fixedFlows.assign(M + 1, 0.0);
    opt.floatWeights.assign(M + 1, 0.0);
    opt.fixedFlows[M] = +accreted(p.rateTime(M));
    for (const Size a : callIdx)
        opt.rights.push_back({a, a, -accreted(p.rateTime(a))});
    FmmLsmConfig cfgB = cfg;
    cfgB.trainingSeed = 333;
    cfgB.valuationSeed = 444444;
    FmmLsmPricer optPricer(model, opt, cfgB);
    const auto optRes = optPricer.calculate();

    const Real lhs = noteRes.lowerBoundCv;
    const Real rhs = noteRes.underlyingValue + optRes.lowerBound;
    const Real se =
        std::sqrt(noteRes.lowerBoundCvSe * noteRes.lowerBoundCvSe + optRes.lowerBoundSe * optRes.lowerBoundSe);
    std::ostringstream probs;
    for (Size r = 0; r < noteRes.exerciseProbability.size(); ++r)
        probs << " " << noteRes.exerciseProbability[r];
    BOOST_TEST_MESSAGE("30NC5 accreting note (cv) " << lhs << " +/- " << noteRes.lowerBoundCvSe << "; zero liability "
                                                    << noteRes.underlyingValue << " + bermudan zero-bond option "
                                                    << optRes.lowerBound << " +/- " << optRes.lowerBoundSe << " = " << rhs
                                                    << "; diff " << lhs - rhs << " vs 3 s.e. " << 3.0 * se << "; runtime "
                                                    << noteRes.runtimeSeconds + optRes.runtimeSeconds << " s");
    BOOST_TEST_MESSAGE("30NC5 call probabilities (5y..29y):" << probs.str() << "; never " << noteRes.noExerciseProbability
                                                            << "; E[notice | call] " << noteRes.expectedExerciseTime);
    BOOST_CHECK_MESSAGE(std::fabs(lhs - rhs) < 3.0 * se, "30NC5 accreting note identity violated: " << lhs - rhs);
    pairedIdentityCheck(notePricer, optPricer, noteRes.underlyingValue, "30NC5 accreting note", optRes.lowerBound);
    {
        const auto dual = notePricer.dualBound(256, 32, 20260925);
        const Real upper = noteRes.lowerBoundCv + dual.gap;
        const Real upperSe = std::sqrt(noteRes.lowerBoundCvSe * noteRes.lowerBoundCvSe + dual.gapSe * dual.gapSe);
        BOOST_TEST_MESSAGE("DUALITY-GAP ROW | callable accreting zero note 30NC5 (issuer) | lower(16k, cv) " << lhs
                           << " +/- " << noteRes.lowerBoundCvSe << " | gap " << dual.gap << " +/- " << dual.gapSe
                           << " | upper = lower + gap " << upper << " +/- " << upperSe << " | outer 256 x inner 32 | "
                           << dual.runtimeSeconds << " s");
        BOOST_CHECK_MESSAGE(dual.gap > -3.0 * dual.gapSe, "30NC5: negative duality gap beyond noise: " << dual.gap);
    }

    // matching accreting cancellable swap: notional accretes, annual fixed coupon accRate on the
    // accreted notional, semiannual float on the accreted notional; cancel rights as the calls
    FmmCallableInstrument accSwap;
    accSwap.style = FmmCallableInstrument::Style::Cancel;
    accSwap.lastFlowIdx = M;
    accSwap.fixedFlows.assign(M + 1, 0.0);
    accSwap.floatWeights.assign(M + 1, 0.0);
    for (Size j = 1; j <= M; ++j)
        accSwap.floatWeights[j] = accreted(p.rateTime(j - 1));
    for (Size c = 2; c <= M; c += 2)
        accSwap.fixedFlows[c] = -accRate * (p.rateTime(c) - p.rateTime(c - 2)) * accreted(p.rateTime(c - 2));
    for (const Size a : callIdx)
        accSwap.rights.push_back({a, a, 0.0});
    FmmLsmPricer accSwapPricer(model, accSwap, cfg);
    const auto accSwapRes = accSwapPricer.calculate();
    FmmCallableInstrument accRecv = accSwap;
    accRecv.style = FmmCallableInstrument::Style::Enter;
    for (Size j = 1; j <= M; ++j) {
        accRecv.floatWeights[j] *= -1.0;
        accRecv.fixedFlows[j] *= -1.0;
    }
    FmmLsmPricer accRecvPricer(model, accRecv, cfgB);
    const auto accRecvRes = accRecvPricer.calculate();
    const Real diffAcc = accSwapRes.lowerBoundCv - (accSwapRes.underlyingValue + accRecvRes.lowerBound);
    const Real seAcc = std::sqrt(accSwapRes.lowerBoundCvSe * accSwapRes.lowerBoundCvSe +
                                 accRecvRes.lowerBoundSe * accRecvRes.lowerBoundSe);
    BOOST_TEST_MESSAGE("30NC5 accreting cancellable swap (cv) " << accSwapRes.lowerBoundCv << " +/- "
                                                                << accSwapRes.lowerBoundCvSe << " vs swap+option "
                                                                << accSwapRes.underlyingValue + accRecvRes.lowerBound
                                                                << "; diff " << diffAcc << " vs 3 s.e. " << 3.0 * seAcc);
    BOOST_CHECK_MESSAGE(std::fabs(diffAcc) < 3.0 * seAcc, "30NC5 accreting swap parity violated: " << diffAcc);
    pairedIdentityCheck(accSwapPricer, accRecvPricer, accSwapRes.underlyingValue, "30NC5 accreting cancellable swap",
                        accRecvRes.lowerBound);
}

BOOST_AUTO_TEST_CASE(testShiftAdmissibilityGuard) {
    BOOST_TEST_MESSAGE("Shifts above 1/tau are rejected at construction (FMM_SPEC.md section 2.2)...");
    Settings::instance().evaluationDate() = Date(19, September, 2026);
    Handle<YieldTermStructure> curve(
        QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.03, Actual365Fixed()));
    Array rateTimes(3);
    rateTimes[0] = 0.0;
    rateTimes[1] = 0.25;
    rateTimes[2] = 0.5;
    Array shifts(2, 4.5); // > 1/tau = 4
    Array volTimes;
    std::vector<Array> volLevels(2, Array(1, 0.01));
    BOOST_CHECK_THROW(FmmParametrization(EURCurrency(), curve, rateTimes, shifts, volTimes, volLevels,
                                         FmmParametrization::LocalVolType::DisplacedDiffusion, 0.6, 0.08, 1),
                      QuantLib::Error);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
