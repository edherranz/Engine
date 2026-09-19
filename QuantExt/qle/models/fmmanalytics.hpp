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

/*! \file fmmanalytics.hpp
    \brief FMM closed forms: backward/forward caplets under displaced-diffusion or normal
           dynamics, and the frozen-gradient normal swaption approximation
           (ore-fmm docs/FMM_SPEC.md sections 7.2-7.4); initial domain per spec section 7.4:
           physically settled, forward-starting, single-curve swaps on tenor-grid schedules
    \ingroup models
*/

#ifndef quantext_fmm_analytics_hpp
#define quantext_fmm_analytics_hpp

#include <qle/models/fmmparametrization.hpp>

#include <ql/option.hpp>

namespace QuantExt {
using namespace QuantLib;

//! on-grid swap: float leg pays tau_j R_j at T_j for j in (a, b]; fixed leg pays
//! fixedAccruals[c] * K at T_{fixedPayIndices[c]}, indices strictly increasing within (a, b]
struct FmmSwapSpec {
    Size a = 0, b = 0;
    std::vector<Size> fixedPayIndices;
    std::vector<Real> fixedAccruals;
    void validate(const Size M) const;
};

struct FmmSwaptionApproxResult {
    Real price = 0.0;      //!< option value (unit notional)
    Real normalVol = 0.0;  //!< implied normal vol of the forward swap rate, annualised over T_a
    Real forward = 0.0;    //!< forward swap rate S(0)
    Real annuity = 0.0;    //!< fixed-leg annuity A(0)
    Real variance = 0.0;   //!< total normal variance of S at expiry
};

//! undiscounted-forward caplet building blocks
Real fmmCapletStdDev(const FmmParametrization& p, const Size j, const bool backwardLooking);

//! caplet/floorlet on period j, unit notional, payoff tau_j (w (R - K))^+ paid at T_j
Real fmmCapletPrice(const FmmParametrization& p, const Size j, const Real K, const bool backwardLooking,
                    const Option::Type type = Option::Call);

//! implied normal vol of the caplet, annualised over the expiry (T_j backward, T_{j-1} forward)
Real fmmCapletNormalVol(const FmmParametrization& p, const Size j, const Real K, const bool backwardLooking);

Real fmmForwardSwapRate(const FmmParametrization& p, const FmmSwapSpec& swap);
Real fmmAnnuity(const FmmParametrization& p, const FmmSwapSpec& swap);

//! frozen-gradient normal approximation (FMM_SPEC.md section 7.4): dS = sum_i q_i dR_i with the
//! exact gradient q_i = dS/dR_i at t=0 (weights and weight-derivatives combined), local vols
//! frozen at time-0 rates, time integrals exact; expiry = T_a
FmmSwaptionApproxResult fmmSwaptionApprox(const FmmParametrization& p, const FmmSwapSpec& swap, const Real K,
                                          const Option::Type type = Option::Call);

} // namespace QuantExt

#endif
