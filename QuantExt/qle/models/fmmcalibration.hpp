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

/*! \file fmmcalibration.hpp
    \brief FMM calibration in the separable form lambda_j(t) = Lambda_j a(t) with piecewise-
           constant a and the identifiability normalization a(0) = 1 (ore-fmm docs/FMM_SPEC.md
           section 6.1): per-bucket caplet level bootstrap, sequential coterminal-swaption
           time-dependence bootstrap, and the joint iteration. The OREData builder/report plumbing
           (model-neutral keys + modelType=FMM + fmm* fields) is the A5 milestone; this layer
           produces the numbers.
    \ingroup models
*/

#ifndef quantext_fmm_calibration_hpp
#define quantext_fmm_calibration_hpp

#include <qle/models/fmmanalytics.hpp>

#include <string>
#include <vector>

namespace QuantExt {
using namespace QuantLib;

struct FmmCapletVolTarget {
    Size bucket = 0;              //!< rate index j
    Real normalVol = 0.0;         //!< target implied normal vol (annualised per fmmCapletNormalVol)
    bool backwardLooking = true;
};

struct FmmSwaptionVolTarget {
    FmmSwapSpec swap;
    Real normalVol = 0.0;
    std::string label;
};

//! separable volatility state; segmentTimes must equal the parametrization's vol breakpoints
struct FmmSeparableVols {
    Array segmentTimes;        //!< breakpoints of a(t)
    std::vector<Real> levels;  //!< Lambda_j, size M
    std::vector<Real> a;       //!< segment values, size segmentTimes.size() + 1
    //! push Lambda_j * a_k into the parametrization's raw parameter storage
    void apply(FmmParametrization& p) const;
    //! rescale to the identifiability normalization a(0) = 1, keeping the model invariant
    void normalize();
};

struct FmmCalibrationReport {
    struct Row {
        std::string instrument;
        Real marketVol = 0.0, modelVol = 0.0; //!< normal vols
        Real errorBp = 0.0;
    };
    std::vector<Row> rows;
    Size iterations = 0;
    Real runtimeSeconds = 0.0;
    bool converged = false;
    Real maxAbsErrorBp() const;
};

//! strategy (a): solve the per-bucket levels to the caplet targets, holding a(t) fixed; buckets
//! without a target inherit the nearest preceding calibrated level (flat extrapolation)
void fmmCapletLevelBootstrap(FmmParametrization& p, FmmSeparableVols& v,
                             const std::vector<FmmCapletVolTarget>& targets);

//! strategy (b): sequential bootstrap of the a(t) segments to coterminal swaption targets
//! (ascending expiries; segmentTimes must contain each target expiry except the last), holding
//! the levels fixed; normalized to a(0) = 1 on exit
void fmmSwaptionTimeDependenceBootstrap(FmmParametrization& p, FmmSeparableVols& v,
                                        const std::vector<FmmSwaptionVolTarget>& targets);

//! strategy (c): joint iteration — caplets fix levels, swaptions fix the time dependence;
//! returns the calibration report with per-instrument residuals, iterations and runtime
FmmCalibrationReport fmmJointBootstrap(FmmParametrization& p, FmmSeparableVols& v,
                                       const std::vector<FmmCapletVolTarget>& capletTargets,
                                       const std::vector<FmmSwaptionVolTarget>& swaptionTargets,
                                       const Size maxIterations = 10, const Real tolBp = 1e-3);

} // namespace QuantExt

#endif
