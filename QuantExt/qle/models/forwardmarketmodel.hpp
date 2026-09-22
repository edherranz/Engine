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

/*! \file forwardmarketmodel.hpp
    \brief generalized Forward Market Model (Lyashenko-Mercurio) core: state evolution under the
           spot or an extended terminal-forward measure, bank account and discount bonds on and
           off the tenor grid via the FMM-fitted Markovian HJM completion; see docs/FMM_SPEC.md
           in the ore-fmm repository
    \ingroup models
*/

#ifndef quantext_forward_market_model_hpp
#define quantext_forward_market_model_hpp

#include <qle/methods/multipathgeneratorbase.hpp>
#include <qle/methods/multipathvariategenerator.hpp>
#include <qle/models/fmmparametrization.hpp>

namespace QuantExt {
using namespace QuantLib;

//! Forward Market Model
/*! Joint evolution of the backward-looking forward rates R_1..R_M on the tenor grid, with

    - drifts under the spot measure Q (bank-account numeraire) or an extended T_k-forward measure
      (FMM_SPEC.md section 2.3), assembled from exact per-step integrated covariances so the
      decay never enters through point sampling;
    - log-Euler stepping in the displaced variable ln(R_j + delta_j) (or plain Euler for the
      normal type) with predictor-corrector treatment of the state-dependent drift factors;
    - the off-grid completion states (the pathwise variances Y_kk and the live period's
      front-stub pair (x, y), FMM_SPEC.md section 4.3) so that B(t) and P(t,T) are available
      for arbitrary t and T.

    The tenor grid must start at T_0 = 0. Simulation step boundaries must include the tenor dates
    they cross (the path generator below enforces this); steps therefore never span a grid date's
    interior, which keeps the completion snapshots exact.

    \ingroup models
*/
class ForwardMarketModel {
public:
    enum class Measure { Spot, TerminalForward };

    //! full per-path model state
    struct State {
        Time t = 0.0;
        Array R;        //!< rates R_1..R_M (frozen automatically after their T_j)
        Array Ykk;      //!< pathwise variances int_0^t varsigma_k^2, k = 1..M (LM2020 eq. 31)
        Real x = 0.0;   //!< live period front-stub state, resets at each grid date
        Array RSnap;    //!< rates observed at the current period start T_{period-1}
        Array YkkSnap;  //!< Y_kk observed at the current period start
        Size period = 1; //!< eta(t): the accruing period index
    };

    //! precomputed, state-independent data for one step [s, e]
    struct StepData {
        Time s = 0.0, e = 0.0;
        Matrix cov;       //!< M x M exact integrated covariances (both decays, correlation)
        Matrix shockSqrt; //!< (M+1) x (M+1) pseudo square root; row M is the raw front-stub shock
        Array levelVar;   //!< per rate: int lambda_k^2 du clipped at T_k (Y_kk increments)
        Real gInt = 0.0;  //!< int g_live du for the front-stub drift
        Size live = 0;    //!< period index the step lies in (0 if beyond the grid)
    };

    explicit ForwardMarketModel(const QuantLib::ext::shared_ptr<FmmParametrization>& parametrization,
                                const Measure measure = Measure::Spot, const Size terminalIndex = Null<Size>());

    const QuantLib::ext::shared_ptr<FmmParametrization>& parametrization() const { return p_; }
    Measure measure() const { return measure_; }
    Size terminalIndex() const { return terminalIndex_; }

    //! state at t = 0 with curve-implied rates
    State initialState() const;
    //! precompute step data for [s, e]; the step must not cross a tenor date's interior
    StepData makeStep(const Time s, const Time e) const;
    //! evolve the state over one step given M+1 iid standard normal draws
    void evolve(State& state, const StepData& step, const Array& z) const;
    //! evolve given the M+1 correlated shocks directly (v = shockSqrt * z). Because integrated
    //! covariances are additive over subintervals, the sum of fine-step shock vectors over a
    //! coarse window is a distributionally exact coarse shock vector — the seam used by the
    //! paired common-random-numbers bias study (FMM_SPEC.md section 9.5/9.7)
    void evolveWithCorrelatedShocks(State& state, const StepData& step, const Array& v) const;
    //! same, with an additional drift per rate over the step (on ln(R_j + delta_j) for the
    //! displaced type, on R_j for the normal type): the quanto adjustment of a foreign model
    //! simulated under another currency's spot measure (A9.5, FMM_SPEC.md section 8)
    void evolveWithCorrelatedShocks(State& state, const StepData& step, const Array& v,
                                    const Array& extraDrift) const;

    //! bank account B(t) (LM2020 eq. 26)
    Real bankAccount(const State& state) const;
    //! extended discount bond P(t, T), T >= t (grid identities + completion stubs)
    Real discountBond(const State& state, const Time T) const;
    //! numeraire: B(t) under Spot, P(t, T_terminal) under TerminalForward
    Real numeraire(const State& state) const;

private:
    //! P(t, T) for t <= T <= T_period via the front stub (LM2020 eq. 32)
    Real frontStub(const State& state, const Time T) const;
    //! forward bond P(t, T_{k-1}, T), t <= T_{k-1} < T <= T_k (LM2020 eq. 30), given R_k, Y_kk at t
    Real backStub(const Size k, const Real Rk, const Real YkkAtT, const Time T) const;
    //! curve as of the period-start snapshot, maturity u in [snap, T_period]
    Real snapCurve(const State& state, const Time u) const;
    Real P0(const Time T) const { return p_->termStructure()->discount(T); }
    //! q_i = tau_i phi_i(R_i) / (1 + tau_i R_i), the frozen drift factors
    void driftFactors(const Array& R, Array& q) const;
    //! drift of ln(R_j+delta) (DD) resp. R_j (Normal) over the step, given factors q
    void addDrift(const StepData& step, const Array& q, Array& drift) const;

    QuantLib::ext::shared_ptr<FmmParametrization> p_;
    Measure measure_;
    Size terminalIndex_;
    Size M_;
};

//! path generator producing full model states on a user time grid
/*! The internal simulation grid is the union of the user times and all tenor dates up to the
    horizon; states are reported at the user times. Uses ORE's variate generators (MersenneTwister,
    Sobol, Sobol + Brownian bridge) with dimension M+1 per step.
    \ingroup models
*/
class FmmPathGenerator {
public:
    struct Path {
        //! states at the requested (user) times, in order
        std::vector<ForwardMarketModel::State> states;
    };

    FmmPathGenerator(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model, const std::vector<Time>& userTimes,
                     const SequenceType sequenceType = SobolBrownianBridge, const BigNatural seed = 42,
                     const SobolBrownianGenerator::Ordering ordering = SobolBrownianGenerator::Steps,
                     const SobolRsg::DirectionIntegers directionIntegers = SobolRsg::JoeKuoD7);

    const std::vector<Time>& userTimes() const { return userTimes_; }
    const std::vector<Time>& simulationTimes() const { return simTimes_; }
    Path next() const;
    void reset();

private:
    QuantLib::ext::shared_ptr<ForwardMarketModel> model_;
    std::vector<Time> userTimes_, simTimes_;
    std::vector<Size> userIndex_; // position of each user time in simTimes_
    std::vector<ForwardMarketModel::StepData> steps_;
    SequenceType sequenceType_;
    BigNatural seed_;
    SobolBrownianGenerator::Ordering ordering_;
    SobolRsg::DirectionIntegers directionIntegers_;
    QuantLib::ext::shared_ptr<MultiPathVariateGeneratorBase> variates_;
};

} // namespace QuantExt

#endif
