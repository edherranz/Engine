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

#include <qle/models/fmmcalibration.hpp>
#include <qle/models/forwardmarketmodel.hpp>

#include <ql/math/comparison.hpp>
#include <ql/math/solvers1d/brent.hpp>
#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/costfunction.hpp>
#include <ql/math/optimization/endcriteria.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/math/optimization/problem.hpp>

#include <chrono>
#include <sstream>

namespace QuantExt {

namespace {
Real capletForward(const FmmParametrization& p, const FmmCapletVolTarget& t) {
    return (p.termStructure()->discount(p.rateTime(t.bucket - 1)) / p.termStructure()->discount(p.rateTime(t.bucket)) -
            1.0) /
           p.tau(t.bucket);
}
Real modelSwaptionVol(const FmmParametrization& p, const FmmSwaptionVolTarget& t,
                      const FmmSwaptionApproxMethod method) {
    return fmmSwaptionApprox(p, t.swap, fmmTargetStrike(p, t), Option::Call, method).normalVol;
}
} // namespace

Real fmmTargetStrike(const FmmParametrization& p, const FmmSwaptionVolTarget& t) {
    return t.strike == Null<Real>() ? fmmForwardSwapRate(p, t.swap) : t.strike;
}

void FmmSeparableVols::apply(FmmParametrization& p) const {
    QL_REQUIRE(levels.size() == p.numberOfRates(), "FmmSeparableVols: levels size mismatch");
    QL_REQUIRE(segmentTimes.size() + 1 == a.size(), "FmmSeparableVols: a size mismatch");
    const Array& pt = p.parameterTimes(0);
    QL_REQUIRE(pt.size() == segmentTimes.size(), "FmmSeparableVols: segmentTimes must equal the "
                                                 "parametrization's vol breakpoints");
    for (Size k = 0; k < pt.size(); ++k)
        QL_REQUIRE(QuantLib::close_enough(pt[k], segmentTimes[k]),
                   "FmmSeparableVols: segment time " << segmentTimes[k] << " != breakpoint " << pt[k]);
    const auto param = p.parameter(0);
    const Size nVals = a.size();
    Size idx = 0;
    for (Size j = 0; j < levels.size(); ++j)
        for (Size k = 0; k < nVals; ++k)
            param->setParam(idx++, levels[j] * a[k]);
    p.update();
}

void FmmSeparableVols::normalize() {
    const Real a0 = a.front();
    QL_REQUIRE(a0 > 0.0, "FmmSeparableVols: non-positive leading segment");
    for (auto& ak : a)
        ak /= a0;
    for (auto& l : levels)
        l *= a0;
}

Real FmmCalibrationReport::maxAbsErrorBp() const {
    Real m = 0.0;
    for (const auto& r : rows)
        m = std::max(m, std::fabs(r.errorBp));
    return m;
}

void fmmCapletLevelBootstrap(FmmParametrization& p, FmmSeparableVols& v,
                             const std::vector<FmmCapletVolTarget>& targets) {
    QL_REQUIRE(!targets.empty(), "fmmCapletLevelBootstrap: no targets");
    std::vector<bool> covered(p.numberOfRates(), false);
    for (const auto& t : targets) {
        QL_REQUIRE(t.bucket >= 1 && t.bucket <= p.numberOfRates(),
                   "fmmCapletLevelBootstrap: bucket " << t.bucket << " out of range");
        QL_REQUIRE(!covered[t.bucket - 1], "fmmCapletLevelBootstrap: duplicate bucket " << t.bucket);
        covered[t.bucket - 1] = true;
        const Real F0 = capletForward(p, t);
        Brent solver;
        v.levels[t.bucket - 1] = solver.solve(
            [&](const Real level) {
                v.levels[t.bucket - 1] = level;
                v.apply(p);
                return fmmCapletNormalVol(p, t.bucket, F0, t.backwardLooking) - t.normalVol;
            },
            1e-12, std::max(v.levels[t.bucket - 1], 1e-4), 1e-8, 50.0);
    }
    // flat extrapolation for uncovered buckets from the nearest preceding calibrated level
    Real last = Null<Real>();
    for (Size j = 0; j < p.numberOfRates(); ++j) {
        if (covered[j])
            last = v.levels[j];
        else if (last != Null<Real>())
            v.levels[j] = last;
    }
    // leading uncovered buckets take the first calibrated level
    for (Size j = p.numberOfRates(); j-- > 0;) {
        if (covered[j])
            last = v.levels[j];
        else if (last != Null<Real>())
            v.levels[j] = last;
    }
    v.apply(p);
}

void fmmSwaptionTimeDependenceBootstrap(FmmParametrization& p, FmmSeparableVols& v,
                                        const std::vector<FmmSwaptionVolTarget>& targets,
                                        const FmmSwaptionApproxMethod method) {
    QL_REQUIRE(!targets.empty(), "fmmSwaptionTimeDependenceBootstrap: no targets");
    for (Size k = 1; k < targets.size(); ++k)
        QL_REQUIRE(targets[k].swap.a > targets[k - 1].swap.a,
                   "fmmSwaptionTimeDependenceBootstrap: expiries must be ascending");
    QL_REQUIRE(targets.size() == v.a.size(),
               "fmmSwaptionTimeDependenceBootstrap: need one segment per target, got "
                   << v.a.size() << " segments for " << targets.size() << " targets");
    for (Size k = 0; k + 1 < targets.size(); ++k) {
        const Time Te = p.rateTime(targets[k].swap.a);
        QL_REQUIRE(QuantLib::close_enough(v.segmentTimes[k], Te),
                   "fmmSwaptionTimeDependenceBootstrap: segment boundary " << v.segmentTimes[k]
                                                                           << " != target expiry " << Te);
    }
    for (Size k = 0; k < targets.size(); ++k) {
        Brent solver;
        v.a[k] = solver.solve(
            [&](const Real ak) {
                v.a[k] = ak;
                v.apply(p);
                return modelSwaptionVol(p, targets[k], method) - targets[k].normalVol;
            },
            1e-12, std::max(v.a[k], 0.1), 1e-8, 100.0);
    }
    v.normalize();
    v.apply(p);
}

namespace {
//! least-squares residuals (bp of normal vol) of the swaption basket in the a(t) segments (log scale)
class TimeDependenceCost : public CostFunction {
public:
    TimeDependenceCost(FmmParametrization& p, FmmSeparableVols& v, const std::vector<FmmSwaptionVolTarget>& targets,
                       const FmmSwaptionApproxMethod method)
        : p_(p), v_(v), targets_(targets), method_(method) {}
    Array values(const Array& x) const override {
        for (Size k = 0; k < v_.a.size(); ++k)
            v_.a[k] = std::exp(x[k]);
        v_.apply(p_);
        Array r(targets_.size());
        for (Size i = 0; i < targets_.size(); ++i)
            r[i] = (modelSwaptionVol(p_, targets_[i], method_) - targets_[i].normalVol) * 1e4;
        return r;
    }
    Real value(const Array& x) const override {
        const Array r = values(x);
        return std::sqrt(DotProduct(r, r) / static_cast<Real>(std::max<Size>(r.size(), 1)));
    }

private:
    FmmParametrization& p_;
    FmmSeparableVols& v_;
    const std::vector<FmmSwaptionVolTarget>& targets_;
    FmmSwaptionApproxMethod method_;
};
} // namespace

Real fmmSwaptionTimeDependenceBestFit(FmmParametrization& p, FmmSeparableVols& v,
                                      const std::vector<FmmSwaptionVolTarget>& targets,
                                      const FmmSwaptionApproxMethod method) {
    QL_REQUIRE(!targets.empty(), "fmmSwaptionTimeDependenceBestFit: no targets");
    QL_REQUIRE(v.a.size() == v.segmentTimes.size() + 1, "fmmSwaptionTimeDependenceBestFit: a / segmentTimes mismatch");
    Array x(v.a.size());
    for (Size k = 0; k < v.a.size(); ++k)
        x[k] = std::log(std::max(v.a[k], 1e-6));
    TimeDependenceCost cost(p, v, targets, method);
    NoConstraint constraint;
    Problem problem(cost, constraint, x);
    LevenbergMarquardt lm(1e-10, 1e-10, 1e-10);
    EndCriteria endCriteria(1000, 100, 1e-12, 1e-12, 1e-12);
    lm.minimize(problem, endCriteria);
    const Array& xs = problem.currentValue();
    for (Size k = 0; k < v.a.size(); ++k)
        v.a[k] = std::exp(xs[k]);
    v.normalize();
    v.apply(p);
    Real ss = 0.0;
    for (const auto& t : targets) {
        const Real r = (modelSwaptionVol(p, t, method) - t.normalVol) * 1e4;
        ss += r * r;
    }
    return std::sqrt(ss / static_cast<Real>(targets.size()));
}

namespace {
FmmCalibrationReport jointAlternation(FmmParametrization& p, FmmSeparableVols& v,
                                      const std::vector<FmmCapletVolTarget>& capletTargets,
                                      const std::vector<FmmSwaptionVolTarget>& swaptionTargets,
                                      const Size maxIterations, const Real tolBp, const FmmSwaptionApproxMethod method,
                                      const bool bestFit) {
    const auto start = std::chrono::steady_clock::now();
    FmmCalibrationReport report;
    // the alternation (caplet levels, then swaption time dependence) is iterated until either all
    // targets are within the tolerance (converged) or the parameters stop moving (stationary: the
    // fixed point of the alternation, reached when the two baskets are not jointly attainable and
    // an irreducible residual remains); stopping at a fixed iteration count instead would make the
    // calibrated state, and every sensitivity through it, depend on where the count cuts the path
    for (Size it = 0; it < maxIterations; ++it) {
        report.iterations = it + 1;
        const std::vector<Real> prevLevels = v.levels, prevA = v.a;
        fmmCapletLevelBootstrap(p, v, capletTargets);
        if (bestFit)
            fmmSwaptionTimeDependenceBestFit(p, v, swaptionTargets, method);
        else
            fmmSwaptionTimeDependenceBootstrap(p, v, swaptionTargets, method);
        Real change = 0.0;
        for (Size j = 0; j < v.levels.size(); ++j)
            change = std::max(change, std::fabs(v.levels[j] - prevLevels[j]) / std::max(std::fabs(prevLevels[j]), 1e-12));
        for (Size k = 0; k < v.a.size(); ++k)
            change = std::max(change, std::fabs(v.a[k] - prevA[k]) / std::max(std::fabs(prevA[k]), 1e-12));
        report.lastParameterChange = change;
        // convergence: all targets repriced within tolerance
        Real worst = 0.0;
        for (const auto& t : capletTargets)
            worst = std::max(worst, std::fabs(fmmCapletNormalVol(p, t.bucket, capletForward(p, t), t.backwardLooking) -
                                              t.normalVol) *
                                        1e4);
        for (const auto& t : swaptionTargets)
            worst = std::max(worst, std::fabs(modelSwaptionVol(p, t, method) - t.normalVol) * 1e4);
        if (worst < tolBp) {
            report.converged = true;
            break;
        }
        // stationarity: the exact bootstrap solves converge to machine precision (1e-9), the
        // least-squares step has a noise floor of about 1e-8 in the parameters from its
        // finite-difference Jacobian, so its fixed point is declared at 1e-7 (the residuals are
        // unchanged to 1e-6 bp at that level)
        if (it > 0 && change < (bestFit ? 1e-7 : 1e-9)) {
            report.stationary = true;
            break;
        }
    }
    for (const auto& t : capletTargets) {
        FmmCalibrationReport::Row row;
        std::ostringstream lbl;
        lbl << "caplet_" << t.bucket << (t.backwardLooking ? "_bwd" : "_fwd");
        row.instrument = lbl.str();
        row.marketVol = t.normalVol;
        row.modelVol = fmmCapletNormalVol(p, t.bucket, capletForward(p, t), t.backwardLooking);
        row.errorBp = (row.modelVol - row.marketVol) * 1e4;
        report.rows.push_back(row);
    }
    for (const auto& t : swaptionTargets) {
        FmmCalibrationReport::Row row;
        row.instrument = t.label.empty() ? "swaption" : t.label;
        row.marketVol = t.normalVol;
        row.modelVol = modelSwaptionVol(p, t, method);
        row.errorBp = (row.modelVol - row.marketVol) * 1e4;
        report.rows.push_back(row);
    }
    report.runtimeSeconds =
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count() *
        1e-6;
    return report;
}
} // namespace

FmmCalibrationReport fmmJointBootstrap(FmmParametrization& p, FmmSeparableVols& v,
                                       const std::vector<FmmCapletVolTarget>& capletTargets,
                                       const std::vector<FmmSwaptionVolTarget>& swaptionTargets,
                                       const Size maxIterations, const Real tolBp,
                                       const FmmSwaptionApproxMethod method) {
    return jointAlternation(p, v, capletTargets, swaptionTargets, maxIterations, tolBp, method, false);
}

FmmCalibrationReport fmmJointBestFit(FmmParametrization& p, FmmSeparableVols& v,
                                     const std::vector<FmmCapletVolTarget>& capletTargets,
                                     const std::vector<FmmSwaptionVolTarget>& swaptionTargets,
                                     const Size maxIterations, const Real tolBp,
                                     const FmmSwaptionApproxMethod method) {
    return jointAlternation(p, v, capletTargets, swaptionTargets, maxIterations, tolBp, method, true);
}

FmmMcCorrectedReport fmmMcCorrectedSwaptionBootstrap(FmmParametrization& p,
                                                     const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                                                     FmmSeparableVols& v,
                                                     const std::vector<FmmSwaptionVolTarget>& targets,
                                                     const FmmSwaptionApproxMethod method, const Size pathsPerRep,
                                                     const Size reps, const BigNatural seed,
                                                     const Size maxIterations, const Real tolBp) {
    const auto start = std::chrono::steady_clock::now();
    QL_REQUIRE(model && model->parametrization().get() == &p,
               "fmmMcCorrectedSwaptionBootstrap: the model must wrap the parametrization being calibrated");
    QL_REQUIRE(!targets.empty(), "fmmMcCorrectedSwaptionBootstrap: no targets");
    FmmMcCorrectedReport report;
    std::vector<FmmSwaptionVolTarget> adjusted = targets;
    for (Size it = 0; it < maxIterations; ++it) {
        report.iterations = it + 1;
        fmmSwaptionTimeDependenceBootstrap(p, v, adjusted, method);
        // same-model Monte Carlo reprice of every market target at its strike (one batch per
        // target: coterminal targets have distinct expiries)
        report.rows.clear();
        Real worst = 0.0;
        for (Size k = 0; k < targets.size(); ++k) {
            const Real K = fmmTargetStrike(p, targets[k]);
            const auto mc = fmmSwaptionMc(model, {targets[k].swap}, {K}, Option::Call, pathsPerRep, reps,
                                          seed + 1000 * static_cast<BigNatural>(k))
                                .front();
            FmmMcCorrectedReport::Row row;
            row.instrument = targets[k].label.empty() ? "swaption" : targets[k].label;
            row.strike = K;
            row.marketVol = targets[k].normalVol;
            row.approxVol = modelSwaptionVol(p, targets[k], method);
            row.mcVol = mc.normalVol;
            row.mcVolSe = mc.normalVolSe;
            row.residualBp = (mc.normalVol - targets[k].normalVol) * 1e4;
            report.rows.push_back(row);
            worst = std::max(worst, std::fabs(row.residualBp));
            // control-variate update of the analytic target by the measured bias
            adjusted[k].normalVol -= (mc.normalVol - targets[k].normalVol);
        }
        report.worstResidualBp.push_back(worst);
        if (worst < tolBp) {
            report.converged = true;
            break;
        }
    }
    report.runtimeSeconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() *
        1e-3;
    return report;
}

} // namespace QuantExt
