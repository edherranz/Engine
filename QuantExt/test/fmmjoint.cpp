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

#include <qle/models/fmmjointmodel.hpp>
#include <qle/models/fmmparametrization.hpp>
#include <qle/models/forwardmarketmodel.hpp>

#include <ql/currencies/america.hpp>
#include <ql/currencies/europe.hpp>
#include <ql/math/statistics/incrementalstatistics.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <cmath>

using namespace QuantLib;
using namespace QuantExt;

BOOST_FIXTURE_TEST_SUITE(QuantExtTestSuite, qle::test::TopLevelFixture)

BOOST_AUTO_TEST_SUITE(FmmJointTest)

namespace {

/* A9.5 test-bed: a domestic (USD, flat 4%) and a foreign (EUR, flat 2.5%) FMM on one quarterly
   five-year grid, displaced diffusion with the 1/tau shift, three factors each, a lognormal FX
   (USD per EUR, spot 1.03) with a piecewise volatility, and the level-factor correlations. */
struct JointBed {
    JointBed(const Real levelD, const Real levelF, const Array& fxVolTimes, const Array& fxVols,
             const FmmJointModel::Correlations& corr, const Size factors = 3) {
        Settings::instance().evaluationDate() = Date(19, September, 2026);
        curveD = Handle<YieldTermStructure>(QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.04, Actual365Fixed()));
        curveF = Handle<YieldTermStructure>(QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), 0.025, Actual365Fixed()));
        const Size M = 20;
        Array rateTimes(M + 1);
        for (Size k = 0; k <= M; ++k)
            rateTimes[k] = 0.25 * static_cast<Real>(k);
        Array shifts(M, 1.0 / 0.25);
        auto make = [&](const Currency& ccy, const Handle<YieldTermStructure>& curve, const Real level) {
            std::vector<Array> volLevels(M, Array(1, level));
            auto p = QuantLib::ext::make_shared<FmmParametrization>(ccy, curve, rateTimes, shifts, Array(), volLevels,
                                                                    FmmParametrization::LocalVolType::DisplacedDiffusion,
                                                                    0.6, 0.08, factors);
            return QuantLib::ext::make_shared<ForwardMarketModel>(p);
        };
        domestic = make(USDCurrency(), curveD, levelD);
        foreign = make(EURCurrency(), curveF, levelF);
        model = QuantLib::ext::make_shared<FmmJointModel>(domestic, foreign, 1.03, fxVolTimes, fxVols, corr);
    }
    Handle<YieldTermStructure> curveD, curveF;
    QuantLib::ext::shared_ptr<ForwardMarketModel> domestic, foreign;
    QuantLib::ext::shared_ptr<FmmJointModel> model;
};

} // namespace

BOOST_AUTO_TEST_CASE(testJointFxMartingaleAndForeignBonds) {
    BOOST_TEST_MESSAGE("A9.5: joint two-currency FMM with FX: X B_f / B_d is a martingale under the domestic spot "
                       "measure at grid and off-grid times, the domestic value of a foreign unit E[X(t)/B_d(t)] equals "
                       "X_0 P_f(0,t) (the quanto drift of the foreign rates), the domestic identities are unchanged, "
                       "and with zero FX vol and zero correlations the FX follows the bank-account ratio pathwise...");
    FmmJointModel::Correlations corr;
    corr.ratesDF = 0.25;
    corr.fxD = -0.2;
    corr.fxF = 0.3;
    JointBed bed(0.0025, 0.0020, Array(1, 2.0), Array({0.08, 0.10}), corr);
    const std::vector<Time> times = {0.25, 0.6, 1.0, 1.3, 2.0, 2.9, 3.0, 4.0};
    const Size paths = 16384;
    FmmJointPathGenerator gen(bed.model, times, SobolBrownianBridge, 4242);
    std::vector<IncrementalStatistics> mart(times.size()), fxDf(times.size()), domDf(times.size());
    for (Size n = 0; n < paths; ++n) {
        const auto path = gen.next();
        for (Size i = 0; i < times.size(); ++i) {
            const auto& st = path.states[i];
            const Real bd = bed.model->domesticBankAccount(st), bf = bed.model->foreignBankAccount(st);
            mart[i].add(st.fx * bf / bd);
            fxDf[i].add(st.fx / bd);
            domDf[i].add(1.0 / bd);
        }
    }
    for (Size i = 0; i < times.size(); ++i) {
        const Real m = mart[i].mean(), se = mart[i].errorEstimate();
        const Real fd = fxDf[i].mean(), fdSe = fxDf[i].errorEstimate();
        const Real target = 1.03 * bed.curveF->discount(times[i]);
        const Real dd = domDf[i].mean(), ddSe = domDf[i].errorEstimate();
        BOOST_TEST_MESSAGE("t = " << times[i] << ": E[X B_f/B_d] = " << m << " +/- " << se << " (X_0 = 1.03), E[X/B_d] = "
                                  << fd << " +/- " << fdSe << " vs X_0 P_f = " << target << " (" << (fd / target - 1.0) * 1e4
                                  << " bp), E[1/B_d] = " << dd << " +/- " << ddSe << " vs P_d = "
                                  << bed.curveD->discount(times[i]));
        BOOST_CHECK_SMALL(m - 1.03, 3.0 * se + 1e-12);
        BOOST_CHECK_SMALL(fd - target, 3.0 * fdSe + 3e-4 * target);
        BOOST_CHECK_SMALL(dd - bed.curveD->discount(times[i]), 3.0 * ddSe + 1e-4);
    }
    // the quanto drift matters: without it the foreign-unit identity fails at these correlations
    FmmJointModel::Correlations none;
    JointBed bed0(0.0025, 0.0020, Array(), Array(1, 0.0), none);
    FmmJointPathGenerator gen0(bed0.model, times, SobolBrownianBridge, 4242);
    Real worst = 0.0;
    for (Size n = 0; n < 512; ++n) {
        const auto path = gen0.next();
        for (Size i = 0; i < times.size(); ++i) {
            const auto& st = path.states[i];
            worst = std::max(worst, std::fabs(st.fx * bed0.model->foreignBankAccount(st) / bed0.model->domesticBankAccount(st) - 1.03));
        }
    }
    BOOST_TEST_MESSAGE("zero FX vol, zero correlations: max |X B_f/B_d - X_0| over paths and times = " << worst);
    BOOST_CHECK_SMALL(worst, 1e-10);
}

BOOST_AUTO_TEST_CASE(testJointShockCorrelations) {
    BOOST_TEST_MESSAGE("A9.5: the joint shock covariance reproduces each model's own covariance block exactly and the "
                       "sample correlations of the log-rate and log-FX increments over one step match the "
                       "level-factor correlations of the inputs...");
    FmmJointModel::Correlations corr;
    corr.ratesDF = 0.25;
    corr.fxD = -0.2;
    corr.fxF = 0.3;
    JointBed bed(0.0025, 0.0020, Array(), Array(1, 0.08), corr);
    const auto step = bed.model->makeStep(0.0, 0.25);
    const Matrix joint = step.shockSqrt * transpose(step.shockSqrt);
    const Matrix covD = step.d.shockSqrt * transpose(step.d.shockSqrt);
    const Matrix covF = step.f.shockSqrt * transpose(step.f.shockSqrt);
    const Size nd = covD.rows(), nf = covF.rows();
    Real worst = 0.0;
    for (Size a = 0; a < nd; ++a)
        for (Size b = 0; b < nd; ++b)
            worst = std::max(worst, std::fabs(joint[a][b] - covD[a][b]));
    for (Size a = 0; a < nf; ++a)
        for (Size b = 0; b < nf; ++b)
            worst = std::max(worst, std::fabs(joint[nd + a][nd + b] - covF[a][b]));
    BOOST_TEST_MESSAGE("max |joint block - own covariance| = " << worst);
    BOOST_CHECK_SMALL(worst, 1e-12);
    // sample correlations over the first step
    const Size paths = 32768, j = 5;
    FmmJointPathGenerator gen(bed.model, {0.25}, SobolBrownianBridge, 99);
    const auto s0 = bed.model->initialState();
    const Real delta = 1.0 / 0.25;
    IncrementalStatistics sD, sF, sX;
    Real cDF = 0.0, cXD = 0.0, cXF = 0.0;
    std::vector<Real> dD, dF, dX;
    for (Size n = 0; n < paths; ++n) {
        const auto st = gen.next().states[0];
        dD.push_back(std::log((st.d.R[j - 1] + delta) / (s0.d.R[j - 1] + delta)));
        dF.push_back(std::log((st.f.R[j - 1] + delta) / (s0.f.R[j - 1] + delta)));
        dX.push_back(std::log(st.fx / s0.fx));
        sD.add(dD.back());
        sF.add(dF.back());
        sX.add(dX.back());
    }
    for (Size n = 0; n < paths; ++n) {
        cDF += (dD[n] - sD.mean()) * (dF[n] - sF.mean());
        cXD += (dX[n] - sX.mean()) * (dD[n] - sD.mean());
        cXF += (dX[n] - sX.mean()) * (dF[n] - sF.mean());
    }
    const Real rDF = cDF / (paths * sD.standardDeviation() * sF.standardDeviation());
    const Real rXD = cXD / (paths * sX.standardDeviation() * sD.standardDeviation());
    const Real rXF = cXF / (paths * sX.standardDeviation() * sF.standardDeviation());
    const Real lD = bed.domestic->parametrization()->factorLoadings()[j - 1][0];
    const Real lF = bed.foreign->parametrization()->factorLoadings()[j - 1][0];
    // expected: the covariance entries of the step (the FX shock and the rates' integrated shocks)
    const Real eDF = joint[j - 1][nd + j - 1] / std::sqrt(joint[j - 1][j - 1] * joint[nd + j - 1][nd + j - 1]);
    const Size x = nd + nf;
    const Real eXD = joint[j - 1][x] / std::sqrt(joint[j - 1][j - 1] * joint[x][x]);
    const Real eXF = joint[nd + j - 1][x] / std::sqrt(joint[nd + j - 1][nd + j - 1] * joint[x][x]);
    BOOST_TEST_MESSAGE("rate 5 d/f: sample " << rDF << " expected " << eDF << " (input " << corr.ratesDF << " x " << lD
                                             << " x " << lF << "); FX/d: " << rXD << " expected " << eXD << "; FX/f: "
                                             << rXF << " expected " << eXF);
    BOOST_CHECK_SMALL(rDF - eDF, 0.02);
    BOOST_CHECK_SMALL(rXD - eXD, 0.02);
    BOOST_CHECK_SMALL(rXF - eXF, 0.02);
    BOOST_CHECK_CLOSE(eDF, corr.ratesDF * lD * lF, 1e-6);
    BOOST_CHECK_CLOSE(eXD, corr.fxD * lD, 1e-6);
    BOOST_CHECK_CLOSE(eXF, corr.fxF * lF, 1e-6);
}

BOOST_AUTO_TEST_CASE(testJointFxVolTermStructure) {
    BOOST_TEST_MESSAGE("A9.5: with vanishing rate volatilities the log-FX variance equals the integrated piecewise FX "
                       "variance, breakpoints included...");
    FmmJointModel::Correlations none;
    JointBed bed(1e-8, 1e-8, Array({1.0, 2.5}), Array({0.08, 0.12, 0.10}), none);
    const std::vector<Time> times = {0.5, 1.0, 1.75, 2.5, 3.5};
    const Size paths = 32768;
    FmmJointPathGenerator gen(bed.model, times, SobolBrownianBridge, 7);
    std::vector<IncrementalStatistics> logs(times.size());
    for (Size n = 0; n < paths; ++n) {
        const auto path = gen.next();
        for (Size i = 0; i < times.size(); ++i) {
            const auto& st = path.states[i];
            // strip the deterministic bank-account ratio (rates are constant here)
            logs[i].add(std::log(st.fx * bed.model->foreignBankAccount(st) / bed.model->domesticBankAccount(st) / 1.03));
        }
    }
    auto integratedVar = [&](const Time T) {
        Real v = 0.0;
        Time prev = 0.0;
        for (const Time k : {1.0, 2.5, T}) {
            const Time hi = std::min(k, T);
            if (hi > prev)
                v += bed.model->fxVol(0.5 * (prev + hi)) * bed.model->fxVol(0.5 * (prev + hi)) * (hi - prev);
            prev = hi;
            if (prev >= T)
                break;
        }
        return v;
    };
    for (Size i = 0; i < times.size(); ++i) {
        const Real var = logs[i].variance(), expected = integratedVar(times[i]);
        BOOST_TEST_MESSAGE("t = " << times[i] << ": log-FX variance " << var << " vs integrated " << expected
                                  << " (mean " << logs[i].mean() << ", expected -var/2 = " << -0.5 * expected << ")");
        BOOST_CHECK_CLOSE(var, expected, 4.0);
        BOOST_CHECK_SMALL(logs[i].mean() + 0.5 * expected, 3.0 * logs[i].errorEstimate() + 1e-4);
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
