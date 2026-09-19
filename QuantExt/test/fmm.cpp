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
#include <ql/math/statistics/incrementalstatistics.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <qle/models/fmmparametrization.hpp>
#include <qle/models/forwardmarketmodel.hpp>
#include <qle/models/irlgm1fconstantparametrization.hpp>
#include <qle/models/lgm.hpp>

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

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
