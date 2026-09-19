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

#include <qle/models/forwardmarketmodel.hpp>

#include <ql/math/comparison.hpp>
#include <ql/math/matrixutilities/pseudosqrt.hpp>

#include <algorithm>
#include <set>

namespace QuantExt {

namespace {
constexpr Real tEps = 1.0e-10;
}

ForwardMarketModel::ForwardMarketModel(const QuantLib::ext::shared_ptr<FmmParametrization>& parametrization,
                                       const Measure measure, const Size terminalIndex)
    : p_(parametrization), measure_(measure), terminalIndex_(terminalIndex), M_(parametrization->numberOfRates()) {
    QL_REQUIRE(p_, "ForwardMarketModel: parametrization is null");
    QL_REQUIRE(QuantLib::close_enough(p_->rateTime(0), 0.0),
               "ForwardMarketModel: tenor grid must start at T_0 = 0, got " << p_->rateTime(0));
    if (measure_ == Measure::TerminalForward) {
        QL_REQUIRE(terminalIndex_ != Null<Size>() && terminalIndex_ >= 1 && terminalIndex_ <= M_,
                   "ForwardMarketModel: terminal index required in [1," << M_ << "]");
    } else {
        terminalIndex_ = Null<Size>();
    }
}

ForwardMarketModel::State ForwardMarketModel::initialState() const {
    State s;
    s.t = 0.0;
    s.R = Array(M_);
    for (Size j = 1; j <= M_; ++j) {
        s.R[j - 1] = (P0(p_->rateTime(j - 1)) / P0(p_->rateTime(j)) - 1.0) / p_->tau(j);
        if (p_->volType() == FmmParametrization::LocalVolType::DisplacedDiffusion) {
            QL_REQUIRE(s.R[j - 1] + p_->shift(j) > 0.0,
                       "ForwardMarketModel: initial rate R_" << j << " = " << s.R[j - 1]
                                                             << " violates R + delta > 0 (delta = "
                                                             << p_->shift(j) << ")");
        }
    }
    s.Ykk = Array(M_, 0.0);
    s.x = 0.0;
    s.RSnap = s.R;
    s.YkkSnap = Array(M_, 0.0);
    s.period = 1;
    return s;
}

ForwardMarketModel::StepData ForwardMarketModel::makeStep(const Time s, const Time e) const {
    QL_REQUIRE(e > s + tEps, "ForwardMarketModel::makeStep: need e > s, got [" << s << "," << e << "]");
    for (Size j = 1; j <= M_; ++j) {
        const Time T = p_->rateTime(j);
        QL_REQUIRE(T <= s + tEps || T >= e - tEps,
                   "ForwardMarketModel::makeStep: step [" << s << "," << e << "] crosses tenor date " << T
                                                          << "; include tenor dates as step boundaries");
    }
    StepData d;
    d.s = s;
    d.e = e;
    const Time mid = 0.5 * (s + e);
    const Size etaMid = p_->eta(mid);
    d.live = etaMid <= M_ ? etaMid : 0;

    d.cov = Matrix(M_, M_, 0.0);
    for (Size i = 1; i <= M_; ++i)
        for (Size j = i; j <= M_; ++j)
            d.cov[i - 1][j - 1] = d.cov[j - 1][i - 1] = p_->integratedCovariance(i, j, s, e);

    d.levelVar = Array(M_, 0.0);
    for (Size k = 1; k <= M_; ++k)
        d.levelVar[k - 1] = p_->integratedLevelVariance(k, s, std::min(e, p_->rateTime(k)));

    Matrix full(M_ + 1, M_ + 1, 0.0);
    for (Size a = 0; a < M_; ++a)
        for (Size b = 0; b < M_; ++b)
            full[a][b] = d.cov[a][b];
    if (d.live > 0) {
        d.gInt = p_->integratedDecay(d.live, s, e);
        for (Size j = 1; j <= M_; ++j) {
            const Real c = p_->integratedCovarianceSingleDecay(d.live, j, s, e);
            full[M_][j - 1] = full[j - 1][M_] = c;
        }
        full[M_][M_] = d.levelVar[d.live - 1];
    }
    d.shockSqrt = pseudoSqrt(full, SalvagingAlgorithm::Spectral);
    return d;
}

void ForwardMarketModel::driftFactors(const Array& R, Array& q) const {
    for (Size i = 1; i <= M_; ++i) {
        const Real den = 1.0 + p_->tau(i) * R[i - 1];
        QL_REQUIRE(den > QL_EPSILON, "ForwardMarketModel: 1 + tau_" << i << " R = " << den
                                                                    << " not positive; rate path exploded");
        q[i - 1] = p_->tau(i) * p_->phi(i, R[i - 1]) / den;
    }
}

void ForwardMarketModel::addDrift(const StepData& step, const Array& q, Array& drift) const {
    const bool dd = p_->volType() == FmmParametrization::LocalVolType::DisplacedDiffusion;
    for (Size j = 1; j <= M_; ++j) {
        Real mu = 0.0;
        if (measure_ == Measure::Spot) {
            // sum from i = 1: dead and not-yet-decaying contributions are handled exactly by the
            // integrated covariances (the eta(t) lower limit of LM2019 eq. 18 is implicit in g)
            for (Size i = 1; i <= j; ++i)
                mu += step.cov[i - 1][j - 1] * q[i - 1];
        } else {
            const Size k = terminalIndex_;
            if (j > k) {
                for (Size i = k + 1; i <= j; ++i)
                    mu += step.cov[i - 1][j - 1] * q[i - 1];
            } else if (j < k) {
                for (Size i = j + 1; i <= k; ++i)
                    mu -= step.cov[i - 1][j - 1] * q[i - 1];
            }
        }
        if (dd)
            mu -= 0.5 * step.cov[j - 1][j - 1]; // Ito term of ln(R_j + delta_j)
        drift[j - 1] += mu;
    }
}

void ForwardMarketModel::evolve(State& state, const StepData& step, const Array& z) const {
    QL_REQUIRE(z.size() == M_ + 1, "ForwardMarketModel::evolve: need " << M_ + 1 << " draws, got " << z.size());
    QL_REQUIRE(std::fabs(state.t - step.s) < 1.0e-8,
               "ForwardMarketModel::evolve: state time " << state.t << " does not match step start " << step.s);
    const bool dd = p_->volType() == FmmParametrization::LocalVolType::DisplacedDiffusion;

    // roll the accruing period at the START of the step: a state resting exactly on a tenor date
    // keeps period = eta(t); the snapshot for the next period is taken here, just before evolving
    while (state.period <= M_ && step.s >= p_->rateTime(state.period) - tEps) {
        state.RSnap = state.R;
        state.YkkSnap = state.Ykk;
        state.x = 0.0;
        ++state.period;
    }
    QL_REQUIRE(state.period <= M_, "ForwardMarketModel::evolve: step starts beyond the tenor grid");

    Array v(M_ + 1, 0.0);
    for (Size a = 0; a <= M_; ++a) {
        Real acc = 0.0;
        for (Size b = 0; b <= M_; ++b)
            acc += step.shockSqrt[a][b] * z[b];
        v[a] = acc;
    }

    // predictor
    Array q(M_), drift1(M_, 0.0);
    driftFactors(state.R, q);
    addDrift(step, q, drift1);
    Array Rt(M_);
    for (Size j = 1; j <= M_; ++j) {
        if (p_->rateTime(j) <= step.s + tEps) {
            Rt[j - 1] = state.R[j - 1]; // frozen after fixing
        } else if (dd) {
            const Real delta = p_->shift(j);
            Rt[j - 1] = (state.R[j - 1] + delta) * std::exp(drift1[j - 1] + v[j - 1]) - delta;
        } else {
            Rt[j - 1] = state.R[j - 1] + drift1[j - 1] + v[j - 1];
        }
    }
    // corrector
    Array q2(M_), drift2(M_, 0.0);
    driftFactors(Rt, q2);
    addDrift(step, q2, drift2);

    const Real yBefore = step.live > 0 ? state.Ykk[step.live - 1] - state.YkkSnap[step.live - 1] : 0.0;
    Real scaleLive = 1.0;
    if (step.live > 0) {
        const Size k = step.live;
        scaleLive = p_->phi(k, state.R[k - 1]) / (state.R[k - 1] + 1.0 / p_->tau(k));
    }

    for (Size j = 1; j <= M_; ++j) {
        if (p_->rateTime(j) <= step.s + tEps)
            continue;
        const Real mu = 0.5 * (drift1[j - 1] + drift2[j - 1]);
        if (dd) {
            const Real delta = p_->shift(j);
            state.R[j - 1] = (state.R[j - 1] + delta) * std::exp(mu + v[j - 1]) - delta;
        } else {
            state.R[j - 1] += mu + v[j - 1];
        }
    }

    // pathwise variance states: Y_kk += int lambda^2 * scale^2 with the local-vol scale
    // scale_k = phi_k(R)/(R + 1/tau_k) frozen at the step start; since (1 + tau R) = tau (R + 1/tau),
    // this equals the drift factor q_k computed from the step-start rates
    for (Size k = 1; k <= M_; ++k) {
        if (step.levelVar[k - 1] <= 0.0)
            continue;
        const Real scale = q[k - 1];
        state.Ykk[k - 1] += step.levelVar[k - 1] * scale * scale;
    }

    if (step.live > 0) {
        const Real yAfter = state.Ykk[step.live - 1] - state.YkkSnap[step.live - 1];
        state.x += step.gInt * 0.5 * (yBefore + yAfter) + scaleLive * v[M_];
    }

    state.t = step.e;
}

Real ForwardMarketModel::backStub(const Size k, const Real Rk, const Real YkkAtT, const Time T) const {
    const Time Ts = p_->rateTime(k - 1), Te = p_->rateTime(k);
    QL_REQUIRE(T >= Ts - tEps && T <= Te + tEps,
               "ForwardMarketModel::backStub: maturity " << T << " outside period " << k);
    const Real gammaT = p_->decay(k, T);      // G_k(T, T_k)
    const Real G = 1.0 - gammaT;              // G_k(T_{k-1}, T)
    const Real fwd0 = P0(Te) / P0(Ts);        // P(0, T_{k-1}, T_k)
    return (P0(T) / P0(Ts)) * std::pow((1.0 + p_->tau(k) * Rk) * fwd0, -G) *
           std::exp(0.5 * G * gammaT * YkkAtT);
}

Real ForwardMarketModel::snapCurve(const State& state, const Time u) const {
    const Size k = state.period;
    const Time snapT = p_->rateTime(k - 1);
    if (u <= snapT + tEps)
        return 1.0;
    return backStub(k, state.RSnap[k - 1], state.YkkSnap[k - 1], u);
}

Real ForwardMarketModel::frontStub(const State& state, const Time T) const {
    const Size k = state.period;
    QL_REQUIRE(k <= M_, "ForwardMarketModel::frontStub: state is beyond the tenor grid");
    QL_REQUIRE(T >= state.t - tEps && T <= p_->rateTime(k) + tEps,
               "ForwardMarketModel::frontStub: maturity " << T << " outside [t, T_period]");
    const Real G = p_->decay(k, state.t) - p_->decay(k, T); // G_k(t, T)
    const Real y = state.Ykk[k - 1] - state.YkkSnap[k - 1];
    return snapCurve(state, T) / snapCurve(state, state.t) * std::exp(-G * state.x - 0.5 * G * G * y);
}

Real ForwardMarketModel::discountBond(const State& state, const Time T) const {
    QL_REQUIRE(T >= state.t - tEps, "ForwardMarketModel::discountBond: T (" << T << ") < t (" << state.t << ")");
    const Size k = state.period;
    QL_REQUIRE(k <= M_, "ForwardMarketModel::discountBond: state is beyond the tenor grid");
    if (T <= p_->rateTime(k) + tEps)
        return frontStub(state, T);
    const Size m = p_->eta(T);
    QL_REQUIRE(m <= M_, "ForwardMarketModel::discountBond: maturity " << T << " beyond the tenor grid");
    Real p = frontStub(state, p_->rateTime(k));
    for (Size j = k + 1; j <= m - 1; ++j)
        p /= 1.0 + p_->tau(j) * state.R[j - 1];
    p *= backStub(m, state.R[m - 1], state.Ykk[m - 1], T);
    return p;
}

Real ForwardMarketModel::bankAccount(const State& state) const {
    const Size k = state.period;
    QL_REQUIRE(k <= M_, "ForwardMarketModel::bankAccount: state is beyond the tenor grid");
    Real b = frontStub(state, p_->rateTime(k));
    for (Size j = 1; j <= k; ++j)
        b *= 1.0 + p_->tau(j) * state.R[j - 1];
    return b;
}

Real ForwardMarketModel::numeraire(const State& state) const {
    if (measure_ == Measure::Spot)
        return bankAccount(state);
    const Time Tk = p_->rateTime(terminalIndex_);
    if (state.t <= Tk + tEps)
        return discountBond(state, Tk);
    // extended bond beyond its maturity: P(t, Tk) = B(t) / B(Tk), with B(Tk) from the fixings
    Real bTk = 1.0;
    for (Size j = 1; j <= terminalIndex_; ++j)
        bTk *= 1.0 + p_->tau(j) * state.R[j - 1];
    return bankAccount(state) / bTk;
}

FmmPathGenerator::FmmPathGenerator(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                                   const std::vector<Time>& userTimes, const SequenceType sequenceType,
                                   const BigNatural seed, const SobolBrownianGenerator::Ordering ordering,
                                   const SobolRsg::DirectionIntegers directionIntegers)
    : model_(model), userTimes_(userTimes), sequenceType_(sequenceType), seed_(seed), ordering_(ordering),
      directionIntegers_(directionIntegers) {
    QL_REQUIRE(model_, "FmmPathGenerator: model is null");
    QL_REQUIRE(!userTimes_.empty(), "FmmPathGenerator: no simulation times");
    for (Size i = 0; i < userTimes_.size(); ++i) {
        QL_REQUIRE(userTimes_[i] > (i == 0 ? 0.0 : userTimes_[i - 1]) + tEps,
                   "FmmPathGenerator: user times must be positive and strictly increasing");
    }
    const auto& p = model_->parametrization();
    std::set<Time> grid(userTimes_.begin(), userTimes_.end());
    const Time horizon = userTimes_.back();
    for (Size j = 1; j <= p->numberOfRates(); ++j)
        if (p->rateTime(j) < horizon - tEps)
            grid.insert(p->rateTime(j));
    simTimes_.assign(grid.begin(), grid.end());
    userIndex_.clear();
    for (const Time t : userTimes_)
        userIndex_.push_back(std::lower_bound(simTimes_.begin(), simTimes_.end(), t - tEps) - simTimes_.begin());
    Time prev = 0.0;
    for (const Time t : simTimes_) {
        steps_.push_back(model_->makeStep(prev, t));
        prev = t;
    }
    reset();
}

void FmmPathGenerator::reset() {
    variates_ = makeMultiPathVariateGenerator(sequenceType_, model_->parametrization()->numberOfRates() + 1,
                                              steps_.size(), seed_, ordering_, directionIntegers_);
}

FmmPathGenerator::Path FmmPathGenerator::next() const {
    const auto sample = variates_->next();
    QL_REQUIRE(sample.value.size() == steps_.size(), "FmmPathGenerator: unexpected variate dimensioning");
    Path path;
    ForwardMarketModel::State state = model_->initialState();
    Size nextUser = 0;
    for (Size i = 0; i < steps_.size(); ++i) {
        model_->evolve(state, steps_[i], sample.value[i]);
        while (nextUser < userIndex_.size() && userIndex_[nextUser] == i) {
            path.states.push_back(state);
            ++nextUser;
        }
    }
    QL_REQUIRE(path.states.size() == userTimes_.size(), "FmmPathGenerator: user time mapping failed");
    return path;
}

} // namespace QuantExt
