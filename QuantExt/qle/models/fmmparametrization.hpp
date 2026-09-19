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

/*! \file fmmparametrization.hpp
    \brief parametrization of the generalized Forward Market Model (Lyashenko-Mercurio),
           see docs/FMM_SPEC.md in the ore-fmm repository for the specification
    \ingroup models
*/

#ifndef quantext_fmm_parametrization_hpp
#define quantext_fmm_parametrization_hpp

#include <qle/models/irlgm1fparametrization.hpp>
#include <qle/models/parametrization.hpp>

#include <ql/math/matrix.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>

namespace QuantExt {
using namespace QuantLib;

//! FMM parametrization
/*! Holds the tenor structure 0 <= T_0 < ... < T_M, the local-volatility specification
    lambda_j(t) * g_j(t) * phi_j(R_j) for each forward rate R_j (j = 1..M), the displaced-diffusion
    shifts delta_j, and the instantaneous correlation with its rank-F factor reduction.

    Two decay modes (FMM_SPEC.md sections 2.1 and 5.1):
    - linear decay g_j(t) = min{ (T_j - t)^+ / (T_j - T_{j-1}), 1 } with piecewise-constant
      per-rate levels lambda_j(t) on a common time grid;
    - LGM-implied decay (replication mode): lambda_j(t) g_j(t) = alpha(t) [H(T_j) - H(max(t,T_{j-1}))],
      shifts delta_j = 1/tau_j, one factor, perfect correlation, which embeds 1F LGM exactly.

    Integrated covariances are computed exactly per smooth piece: closed form for the
    piecewise-constant/linear-decay case, high-order Gauss-Legendre per smooth piece in
    replication mode (exact for piecewise-linear H, machine-precision otherwise).

    \ingroup models
*/
class FmmParametrization : public Parametrization {
public:
    enum class LocalVolType { DisplacedDiffusion, Normal };

    //! general constructor, linear decay
    FmmParametrization(const Currency& currency, const Handle<YieldTermStructure>& termStructure,
                       const Array& rateTimes, const Array& shifts, const Array& volTimes,
                       const std::vector<Array>& volLevels, const LocalVolType volType, const Real rhoInf,
                       const Real beta, const Size factors, const std::string& name = "FMM");

    //! replication-mode constructor: exact embedding of the given 1F LGM parametrization
    FmmParametrization(const Currency& currency, const Handle<YieldTermStructure>& termStructure,
                       const Array& rateTimes, const QuantLib::ext::shared_ptr<IrLgm1fParametrization>& lgm,
                       const std::string& name = "FMM-LGM1F");

    //! \name Inspectors
    //@{
    Size numberOfRates() const { return M_; }
    Size factors() const { return F_; }
    LocalVolType volType() const { return volType_; }
    bool replicationMode() const { return lgm_ != nullptr; }
    const Handle<YieldTermStructure>& termStructure() const { return termStructure_; }
    const QuantLib::ext::shared_ptr<IrLgm1fParametrization>& lgm() const { return lgm_; }
    //! T_j for j = 0..M
    Time rateTime(const Size j) const;
    //! year fraction tau_j of [T_{j-1}, T_j], j = 1..M
    Real tau(const Size j) const;
    //! displaced-diffusion shift delta_j, j = 1..M
    Real shift(const Size j) const;
    //! smallest j in 1..M with T_j >= t; returns M+1 if t > T_M
    Size eta(const Time t) const;
    //! factor loading matrix C (M x F, unit-norm rows)
    const Matrix& factorLoadings() const { return loadings_; }
    //! effective instantaneous correlation (C C^T)_{ij}, 1-based rate indices
    Real correlation(const Size i, const Size j) const;
    //@}

    //! \name Model functions (1-based rate index j in [1, M])
    //@{
    //! level lambda_j(t), excluding decay and local-vol function
    Real volLevel(const Size j, const Time t) const;
    //! decay g_j(t) in [0, 1]
    Real decay(const Size j, const Time t) const;
    //! local vol function phi_j(r) = r + delta_j (DD) or 1 (Normal)
    Real phi(const Size j, const Real r) const;
    //! exact int_s^e lambda_i(u) g_i(u) lambda_j(u) g_j(u) du * correlation(i, j)
    Real integratedCovariance(const Size i, const Size j, const Time s, const Time e) const;
    //! exact int_s^e lambda_i(u) lambda_j(u) g_j(u) du * correlation(i, j)  (decay on j only;
    //! used for the front-stub state x whose HJM factor carries no decay, FMM_SPEC.md section 4.3)
    Real integratedCovarianceSingleDecay(const Size i, const Size j, const Time s, const Time e) const;
    //! exact int_s^e lambda_k(u)^2 du (no decay; the Y_kk integrand for shift 1/tau, LM2020 eq. 31)
    Real integratedLevelVariance(const Size k, const Time s, const Time e) const;
    //! exact int_s^e g_k(u) du (front-stub drift term, FMM_SPEC.md section 4.3)
    Real integratedDecay(const Size k, const Time s, const Time e) const;
    //! points in (s, e) where the integrand changes analytic form (vol breakpoints, grid dates)
    std::vector<Time> integrationKnots(const Time s, const Time e) const;
    //@}

    //! \name Parametrization interface (calibration surface, refined in A3)
    //@{
    Size numberOfParameters() const override { return lgm_ == nullptr ? 1 : 0; }
    const Array& parameterTimes(const Size i) const override;
    const QuantLib::ext::shared_ptr<Parameter> parameter(const Size i) const override;
    void update() const override;
    //@}

private:
    void checkGrid() const;
    void makeFactorLoadings(const Real rhoInf, const Real beta);
    //! covariance over one smooth piece [s, e] (no knots strictly inside), without correlation;
    //! decayOnI / decayOnJ select which legs carry the decay factor
    Real intervalCovariance(const Size i, const Size j, const Time s, const Time e, const bool decayOnI,
                            const bool decayOnJ) const;

    Handle<YieldTermStructure> termStructure_;
    Array rateTimes_;                   // T_0 .. T_M
    Array taus_;                        // tau_1 .. tau_M
    Array shifts_;                      // delta_1 .. delta_M
    Array volTimes_;                    // breakpoints for piecewise-constant levels
    mutable std::vector<Array> volLevels_; // per rate, volTimes_.size() + 1 values
    LocalVolType volType_;
    Size M_, F_;
    Matrix loadings_;                   // M x F
    Matrix effectiveCorrelation_;       // C C^T
    QuantLib::ext::shared_ptr<IrLgm1fParametrization> lgm_; // replication mode if set
    QuantLib::ext::shared_ptr<Parameter> volParameter_;     // concatenated raw levels (general mode)
};

} // namespace QuantExt

#endif
