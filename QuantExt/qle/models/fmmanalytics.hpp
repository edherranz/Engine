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
           dynamics, the frozen-gradient swaption approximations (strike-flat normal and the
           strike-aware effective-shift variant, ore-fmm docs/FMM_SPEC.md sections 7.2-7.4), and
           the same-model Monte Carlo swaption reference with replication-based uncertainty.
           Initial domain per spec section 7.4: physically settled, forward-starting,
           single-curve swaps on tenor-grid schedules.
    \ingroup models
*/

#ifndef quantext_fmm_analytics_hpp
#define quantext_fmm_analytics_hpp

#include <qle/models/fmmparametrization.hpp>

#include <ql/option.hpp>
#include <ql/utilities/null.hpp>

namespace QuantExt {
using namespace QuantLib;

class ForwardMarketModel;

//! on-grid swap: float leg pays tau_j R_j at T_j for j in (a, b] plus, optionally, deterministic
//! basis amounts basisAmounts[c] (per unit notional) at T_{basisPayIndices[c]}: a term-rate
//! floating leg (EURIBOR over ESTR, say) is the RFR par leg plus its time-0 forward basis, frozen
//! (FMM_SPEC.md section 7.4, deterministic basis); fixed leg pays fixedAccruals[c] * K at
//! T_{fixedPayIndices[c]}; all indices strictly increasing within (a, b]
struct FmmSwapSpec {
    Size a = 0, b = 0;
    std::vector<Size> fixedPayIndices;
    std::vector<Real> fixedAccruals;
    std::vector<Size> basisPayIndices;
    std::vector<Real> basisAmounts;
    void validate(const Size M) const;
};

//! frozen-gradient approximation family (FMM_SPEC.md section 7.4)
enum class FmmSwaptionApproxMethod {
    StrikeFlatNormal, //!< local vols frozen at time-0 rates: normal swap rate, strike-independent vol
    EffectiveShift    //!< local vol projected onto the swap rate: displaced-lognormal swap rate
};

struct FmmSwaptionApproxResult {
    Real price = 0.0;      //!< option value (unit notional)
    Real normalVol = 0.0;  //!< implied normal vol of the forward swap rate at the strike, annualised over T_a
    Real forward = 0.0;    //!< forward swap rate S(0)
    Real annuity = 0.0;    //!< fixed-leg annuity A(0)
    Real variance = 0.0;   //!< normal variance of S at expiry consistent with normalVol
    FmmSwaptionApproxMethod method = FmmSwaptionApproxMethod::StrikeFlatNormal;
    //! effective shift delta_S of the displaced-lognormal swap rate (EffectiveShift method);
    //! Null when the normal limit applies (normal dynamics or a vanishing projected slope)
    Real effectiveShift = Null<Real>();
    Real shiftedStdDev = 0.0; //!< total lognormal std dev of S + delta_S over [0, T_a] (EffectiveShift)
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

//! frozen-gradient approximations (FMM_SPEC.md section 7.4): dS = sum_i q_i dR_i with the exact
//! gradient q_i = dS/dR_i at t=0 (weights and weight-derivatives combined), time integrals
//! exact, expiry = T_a. StrikeFlatNormal freezes the local vols at time-0 rates (normal swap
//! rate). EffectiveShift projects the local vol onto the swap rate: with beta_i =
//! Cov(dR_i, dS)/Var(dS) (integrated covariances) the local variance slope dV/dS at S(0) gives a
//! displaced-lognormal S with shift delta_S = V(0)/D - S(0), D = sum_ik q_i q_k beta_i phi_k
//! IC_ik; for a common rate shift and perfectly correlated rates delta_S reduces to that shift,
//! and normal dynamics are the delta_S -> infinity (normal) limit. Strikes with K + delta_S <= 0
//! are outside the domain and rejected.
FmmSwaptionApproxResult fmmSwaptionApprox(const FmmParametrization& p, const FmmSwapSpec& swap, const Real K,
                                          const Option::Type type = Option::Call,
                                          const FmmSwaptionApproxMethod method = FmmSwaptionApproxMethod::StrikeFlatNormal);

//! same-model Monte Carlo reference: physically settled swaptions sharing the expiry T_a
//! (swaps[i] struck at strikes[i]), randomized Sobol with `reps` independent scrambles of
//! `pathsPerRep` paths each; uncertainty from the replication means (normal vol: one-sided
//! implied-vol difference at one price standard error)
struct FmmSwaptionMcResult {
    Real price = 0.0, priceSe = 0.0;
    Real normalVol = 0.0, normalVolSe = 0.0;
};
std::vector<FmmSwaptionMcResult> fmmSwaptionMc(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                                               const std::vector<FmmSwapSpec>& swaps, const std::vector<Real>& strikes,
                                               const Option::Type type, const Size pathsPerRep, const Size reps,
                                               const BigNatural seed);

} // namespace QuantExt

#endif
