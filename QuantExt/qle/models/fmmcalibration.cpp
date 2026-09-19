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

#include <ql/math/comparison.hpp>
#include <ql/math/solvers1d/brent.hpp>

#include <chrono>
#include <sstream>

namespace QuantExt {

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
        const Real F0 = (p.termStructure()->discount(p.rateTime(t.bucket - 1)) /
                             p.termStructure()->discount(p.rateTime(t.bucket)) -
                         1.0) /
                        p.tau(t.bucket);
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
                                        const std::vector<FmmSwaptionVolTarget>& targets) {
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
                return fmmSwaptionApprox(p, targets[k].swap,
                                         fmmForwardSwapRate(p, targets[k].swap))
                           .normalVol -
                       targets[k].normalVol;
            },
            1e-12, std::max(v.a[k], 0.1), 1e-8, 100.0);
    }
    v.normalize();
    v.apply(p);
}

FmmCalibrationReport fmmJointBootstrap(FmmParametrization& p, FmmSeparableVols& v,
                                       const std::vector<FmmCapletVolTarget>& capletTargets,
                                       const std::vector<FmmSwaptionVolTarget>& swaptionTargets,
                                       const Size maxIterations, const Real tolBp) {
    const auto start = std::chrono::steady_clock::now();
    FmmCalibrationReport report;
    for (Size it = 0; it < maxIterations; ++it) {
        report.iterations = it + 1;
        fmmCapletLevelBootstrap(p, v, capletTargets);
        fmmSwaptionTimeDependenceBootstrap(p, v, swaptionTargets);
        // convergence: all targets repriced within tolerance
        Real worst = 0.0;
        for (const auto& t : capletTargets) {
            const Real F0 = (p.termStructure()->discount(p.rateTime(t.bucket - 1)) /
                                 p.termStructure()->discount(p.rateTime(t.bucket)) -
                             1.0) /
                            p.tau(t.bucket);
            worst = std::max(worst,
                             std::fabs(fmmCapletNormalVol(p, t.bucket, F0, t.backwardLooking) - t.normalVol) * 1e4);
        }
        for (const auto& t : swaptionTargets)
            worst = std::max(
                worst,
                std::fabs(fmmSwaptionApprox(p, t.swap, fmmForwardSwapRate(p, t.swap)).normalVol - t.normalVol) * 1e4);
        if (worst < tolBp) {
            report.converged = true;
            break;
        }
    }
    for (const auto& t : capletTargets) {
        const Real F0 = (p.termStructure()->discount(p.rateTime(t.bucket - 1)) /
                             p.termStructure()->discount(p.rateTime(t.bucket)) -
                         1.0) /
                        p.tau(t.bucket);
        FmmCalibrationReport::Row row;
        std::ostringstream lbl;
        lbl << "caplet_" << t.bucket << (t.backwardLooking ? "_bwd" : "_fwd");
        row.instrument = lbl.str();
        row.marketVol = t.normalVol;
        row.modelVol = fmmCapletNormalVol(p, t.bucket, F0, t.backwardLooking);
        row.errorBp = (row.modelVol - row.marketVol) * 1e4;
        report.rows.push_back(row);
    }
    for (const auto& t : swaptionTargets) {
        FmmCalibrationReport::Row row;
        row.instrument = t.label.empty() ? "swaption" : t.label;
        row.marketVol = t.normalVol;
        row.modelVol = fmmSwaptionApprox(p, t.swap, fmmForwardSwapRate(p, t.swap)).normalVol;
        row.errorBp = (row.modelVol - row.marketVol) * 1e4;
        report.rows.push_back(row);
    }
    report.runtimeSeconds =
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count() *
        1e-6;
    return report;
}

} // namespace QuantExt
