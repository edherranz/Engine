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

/*! \file qle/models/fmmjointmodel.hpp
    \brief joint two-currency FMM with a lognormal FX on one tenor grid (A9.5, ore-fmm
           docs/A9_5_JOINT_MODEL_PLAN.md section 2; derivations in FMM_SPEC.md section 8)

    A domestic FMM (currency d) and a foreign FMM (currency f) share the tenor grid; the FX X(t)
    is the number of domestic units per foreign unit. Everything is simulated under the domestic
    spot measure: the domestic model evolves as in the single-currency case; the foreign model's
    rates carry the quanto drift -rho(X, R^f_j) sigma_X lambda_j g_j dt on their driving Brownian
    motions (change of numeraire from B_f X to B_d); the FX follows the exact discrete recursion
    X(e) = X(s) [B_d(e)/B_d(s)] / [B_f(e)/B_f(s)] exp(v_X - sigma_X^2 (e - s) / 2), so that the
    domestic value of the foreign bank account, X B_f / B_d, is a martingale. Correlations: the
    two models' factor structures stay as they are; across blocks the level (first) factors of the
    two currencies carry rho(d, f), and the FX shock carries rho(X, d) and rho(X, f) against the
    level factors; per step the joint shock covariance (both models' rate and front-stub shocks
    plus the FX) is assembled from the parametrizations' exact integrals and cross integrals and
    factorised by a spectral pseudo square root, as the single-currency model does.
*/

#pragma once

#include <qle/models/forwardmarketmodel.hpp>

#include <ql/math/matrix.hpp>

namespace QuantExt {
using namespace QuantLib;

class FmmJointModel {
public:
    //! correlations between the level factors of the two currencies and the FX
    struct Correlations {
        Real ratesDF = 0.0; //!< rho(R^d, R^f) between the level factors
        Real fxD = 0.0;     //!< rho(X, R^d), X in domestic per foreign
        Real fxF = 0.0;     //!< rho(X, R^f)
    };
    struct State {
        ForwardMarketModel::State d, f;
        Real fx = 0.0; //!< X(t), domestic units per foreign unit
        Time t = 0.0;
    };
    struct StepData {
        ForwardMarketModel::StepData d, f;
        Time s = 0.0, e = 0.0;
        Real fxVariance = 0.0; //!< sigma_X^2 (e - s)
        Matrix shockSqrt;      //!< pseudo square root of the joint shock covariance (dimension() square)
        Array quantoDrift;     //!< per foreign rate: the drift on ln(R + delta) (or R) over the step
    };
    /*! fxVolTimes / fxVols: piecewise-constant FX volatility (n breakpoints, n + 1 values);
        fxSpot: X(0) in domestic units per foreign unit */
    FmmJointModel(const QuantLib::ext::shared_ptr<ForwardMarketModel>& domestic,
                  const QuantLib::ext::shared_ptr<ForwardMarketModel>& foreign, const Real fxSpot,
                  const Array& fxVolTimes, const Array& fxVols, const Correlations& correlations);

    const QuantLib::ext::shared_ptr<ForwardMarketModel>& domestic() const { return d_; }
    const QuantLib::ext::shared_ptr<ForwardMarketModel>& foreign() const { return f_; }
    Real fxSpot() const { return fxSpot_; }
    const Correlations& correlations() const { return corr_; }
    Real fxVol(const Time t) const;
    const Array& fxVolTimes() const { return fxVolTimes_; }
    //! number of iid draws per step: (M_d + 1) + (M_f + 1) + 1
    Size dimension() const { return Md_ + Mf_ + 3; }

    State initialState() const;
    //! the step must not cross a tenor date's interior nor an FX volatility breakpoint
    StepData makeStep(const Time s, const Time e) const;
    void evolve(State& state, const StepData& step, const Array& z) const;

    //! X(t)
    Real fx(const State& state) const { return state.fx; }
    Real domesticBankAccount(const State& state) const { return d_->bankAccount(state.d); }
    Real foreignBankAccount(const State& state) const { return f_->bankAccount(state.f); }
    Real domesticDiscountBond(const State& state, const Time T) const { return d_->discountBond(state.d, T); }
    Real foreignDiscountBond(const State& state, const Time T) const { return f_->discountBond(state.f, T); }
    //! numeraire: the domestic bank account
    Real numeraire(const State& state) const { return d_->bankAccount(state.d); }

    //! int_s^e lambda^a_i g^a_i lambda^b_j g^b_j du across two parametrizations (Gauss-Legendre per
    //! smooth piece of the union of both knot sets; exact for the linear-decay form), decays optional
    static Real crossIntegral(const FmmParametrization& a, const Size i, const bool decayI, const FmmParametrization& b,
                              const Size j, const bool decayJ, const Time s, const Time e);
    //! int_s^e lambda^a_i (g^a_i) du (the covariance with a constant-vol factor per unit of its vol)
    static Real volIntegral(const FmmParametrization& a, const Size i, const bool decayI, const Time s, const Time e);

private:
    QuantLib::ext::shared_ptr<ForwardMarketModel> d_, f_;
    Real fxSpot_;
    Array fxVolTimes_, fxVols_;
    Correlations corr_;
    Size Md_, Mf_;
};

//! joint path generator: states at the user times, steps at every tenor date, FX vol breakpoint
//! and user time
class FmmJointPathGenerator {
public:
    struct Path {
        std::vector<FmmJointModel::State> states;
    };
    FmmJointPathGenerator(const QuantLib::ext::shared_ptr<FmmJointModel>& model, const std::vector<Time>& userTimes,
                          const SequenceType sequenceType = SobolBrownianBridge, const BigNatural seed = 42,
                          const SobolBrownianGenerator::Ordering ordering = SobolBrownianGenerator::Steps,
                          const SobolRsg::DirectionIntegers directionIntegers = SobolRsg::JoeKuoD7);
    const std::vector<Time>& userTimes() const { return userTimes_; }
    const std::vector<Time>& simulationTimes() const { return simTimes_; }
    Path next() const;
    void reset();

private:
    QuantLib::ext::shared_ptr<FmmJointModel> model_;
    std::vector<Time> userTimes_, simTimes_;
    std::vector<Size> userIndex_;
    std::vector<FmmJointModel::StepData> steps_;
    SequenceType sequenceType_;
    BigNatural seed_;
    SobolBrownianGenerator::Ordering ordering_;
    SobolRsg::DirectionIntegers directionIntegers_;
    QuantLib::ext::shared_ptr<MultiPathVariateGeneratorBase> variates_;
};

} // namespace QuantExt
