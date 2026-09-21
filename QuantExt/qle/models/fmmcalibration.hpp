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
           time-dependence bootstrap (ATM or deal strikes, strike-flat or effective-shift
           approximation), the joint iteration, and the MC-corrected swaption bootstrap that
           iterates the analytic fit against the same-model Monte Carlo reprice. The OREData
           builder/report plumbing (model-neutral keys + modelType=FMM + fmm* fields) is the A5
           milestone; this layer produces the numbers.
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
    Real normalVol = 0.0;        //!< target implied normal vol at the strike
    Real strike = Null<Real>();  //!< Null = at the money (forward swap rate at pricing time)
    std::string label;
};

//! the strike a target is priced at (its own strike, or the forward swap rate when Null)
Real fmmTargetStrike(const FmmParametrization& p, const FmmSwaptionVolTarget& t);

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
    bool converged = false;  //!< all targets within the tolerance
    bool stationary = false; //!< the alternation reached its fixed point (parameters stopped moving)
    Real lastParameterChange = Null<Real>(); //!< max relative parameter change in the last iteration
    Real maxAbsErrorBp() const;
};

//! strategy (a): solve the per-bucket levels to the caplet targets, holding a(t) fixed; buckets
//! without a target inherit the nearest preceding calibrated level (flat extrapolation)
void fmmCapletLevelBootstrap(FmmParametrization& p, FmmSeparableVols& v,
                             const std::vector<FmmCapletVolTarget>& targets);

//! strategy (b): sequential bootstrap of the a(t) segments to coterminal swaption targets
//! (ascending expiries; segmentTimes must contain each target expiry except the last), holding
//! the levels fixed; each target is matched in normal-vol terms at its own strike with the given
//! approximation; normalized to a(0) = 1 on exit
void fmmSwaptionTimeDependenceBootstrap(
    FmmParametrization& p, FmmSeparableVols& v, const std::vector<FmmSwaptionVolTarget>& targets,
    const FmmSwaptionApproxMethod method = FmmSwaptionApproxMethod::StrikeFlatNormal);

//! strategy (c): joint iteration — caplets fix levels, swaptions fix the time dependence;
//! returns the calibration report with per-instrument residuals, iterations and runtime
FmmCalibrationReport fmmJointBootstrap(FmmParametrization& p, FmmSeparableVols& v,
                                       const std::vector<FmmCapletVolTarget>& capletTargets,
                                       const std::vector<FmmSwaptionVolTarget>& swaptionTargets,
                                       const Size maxIterations = 10, const Real tolBp = 1e-3,
                                       const FmmSwaptionApproxMethod method = FmmSwaptionApproxMethod::StrikeFlatNormal);

//! strategy (d): least-squares fit of the a(t) segments to an arbitrary swaption basket (an
//! expiry x term grid, several terms per expiry, ...), holding the levels fixed: Levenberg-
//! Marquardt on log a_k with residuals in basis points of normal vol at each target's own strike;
//! normalized to a(0) = 1 on exit (the levels absorb the scale, the model is invariant). Segments
//! no target depends on keep their values. Returns the root-mean-square residual in bp.
Real fmmSwaptionTimeDependenceBestFit(FmmParametrization& p, FmmSeparableVols& v,
                                      const std::vector<FmmSwaptionVolTarget>& targets,
                                      const FmmSwaptionApproxMethod method);

//! strategy (c'): the joint alternation of (c) with the least-squares fit (d) in place of the
//! sequential bootstrap: the caplets fix the levels, the swaption basket (any expiry x term set)
//! fixes the time dependence in the least-squares sense; same stopping rules (all targets within
//! tolBp, or stationary parameters) and the same report as (c). With a grid basket the residual
//! criterion is normally not met and the fixed point is the calibrated state.
FmmCalibrationReport fmmJointBestFit(FmmParametrization& p, FmmSeparableVols& v,
                                     const std::vector<FmmCapletVolTarget>& capletTargets,
                                     const std::vector<FmmSwaptionVolTarget>& swaptionTargets,
                                     const Size maxIterations, const Real tolBp,
                                     const FmmSwaptionApproxMethod method);

//! MC-corrected strategy (b): the analytic bootstrap is iterated against targets shifted by the
//! Monte Carlo-measured approximation bias (adjusted_k <- adjusted_k - (mcVol_k - marketVol_k))
//! until the same-model MC reprices every market target within tolBp or maxIterations is hit.
//! The MC reference uses randomized Sobol replications (fmmSwaptionMc); residuals are reported
//! with their replication standard errors, so a residual is only meaningful relative to them.
struct FmmMcCorrectedReport {
    struct Row {
        std::string instrument;
        Real strike = 0.0;
        Real marketVol = 0.0, approxVol = 0.0, mcVol = 0.0, mcVolSe = 0.0; //!< normal vols
        Real residualBp = 0.0;                                            //!< (mcVol - marketVol) in bp
    };
    std::vector<Row> rows;             //!< final iteration
    std::vector<Real> worstResidualBp; //!< per iteration
    Size iterations = 0;
    bool converged = false;
    Real runtimeSeconds = 0.0;
};
FmmMcCorrectedReport fmmMcCorrectedSwaptionBootstrap(FmmParametrization& p,
                                                     const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                                                     FmmSeparableVols& v,
                                                     const std::vector<FmmSwaptionVolTarget>& targets,
                                                     const FmmSwaptionApproxMethod method, const Size pathsPerRep,
                                                     const Size reps, const BigNatural seed,
                                                     const Size maxIterations = 5, const Real tolBp = 0.2);

} // namespace QuantExt

#endif
