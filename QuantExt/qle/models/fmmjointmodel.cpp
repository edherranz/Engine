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

#include <qle/models/fmmjointmodel.hpp>

#include <ql/math/comparison.hpp>
#include <ql/math/integrals/gaussianquadratures.hpp>
#include <ql/math/matrixutilities/pseudosqrt.hpp>

#include <algorithm>
#include <cmath>
#include <set>

namespace QuantExt {

namespace {
constexpr Real tEps = 1.0e-10;
}

FmmJointModel::FmmJointModel(const QuantLib::ext::shared_ptr<ForwardMarketModel>& domestic,
                             const QuantLib::ext::shared_ptr<ForwardMarketModel>& foreign, const Real fxSpot,
                             const Array& fxVolTimes, const Array& fxVols, const Correlations& correlations)
    : d_(domestic), f_(foreign), fxSpot_(fxSpot), fxVolTimes_(fxVolTimes), fxVols_(fxVols), corr_(correlations) {
    QL_REQUIRE(d_ && f_, "FmmJointModel: models are null");
    QL_REQUIRE(d_->measure() == ForwardMarketModel::Measure::Spot && f_->measure() == ForwardMarketModel::Measure::Spot,
               "FmmJointModel: both models must be set up under their spot measures");
    const auto& pd = *d_->parametrization();
    const auto& pf = *f_->parametrization();
    QL_REQUIRE(pd.currency() != pf.currency(), "FmmJointModel: the two models must be in different currencies");
    Md_ = pd.numberOfRates();
    Mf_ = pf.numberOfRates();
    QL_REQUIRE(Md_ == Mf_, "FmmJointModel: the two models must share the tenor grid (" << Md_ << " vs " << Mf_ << " periods)");
    for (Size j = 0; j <= Md_; ++j)
        QL_REQUIRE(std::fabs(pd.rateTime(j) - pf.rateTime(j)) < 1e-10,
                   "FmmJointModel: tenor date " << j << " differs between the models (" << pd.rateTime(j) << " vs "
                                                << pf.rateTime(j) << ")");
    QL_REQUIRE(fxSpot_ > 0.0, "FmmJointModel: FX spot must be positive");
    QL_REQUIRE(fxVols_.size() == fxVolTimes_.size() + 1, "FmmJointModel: " << fxVols_.size() << " FX vols for "
                                                                             << fxVolTimes_.size() << " breakpoints");
    for (Size i = 0; i < fxVolTimes_.size(); ++i)
        QL_REQUIRE(fxVolTimes_[i] > (i == 0 ? 0.0 : fxVolTimes_[i - 1]), "FmmJointModel: FX vol times must be increasing");
    for (const Real v : fxVols_)
        QL_REQUIRE(v >= 0.0, "FmmJointModel: negative FX vol");
    for (const Real r : {corr_.ratesDF, corr_.fxD, corr_.fxF})
        QL_REQUIRE(r >= -1.0 && r <= 1.0, "FmmJointModel: correlation " << r << " outside [-1, 1]");
}

Real FmmJointModel::fxVol(const Time t) const {
    Size k = 0;
    while (k < fxVolTimes_.size() && t >= fxVolTimes_[k] - tEps)
        ++k;
    return fxVols_[k];
}

FmmJointModel::State FmmJointModel::initialState() const {
    State s;
    s.d = d_->initialState();
    s.f = f_->initialState();
    s.fx = fxSpot_;
    s.t = 0.0;
    return s;
}

Real FmmJointModel::crossIntegral(const FmmParametrization& a, const Size i, const bool decayI,
                                  const FmmParametrization& b, const Size j, const bool decayJ, const Time s, Time e) {
    if (decayI)
        e = std::min(e, a.rateTime(i));
    if (decayJ)
        e = std::min(e, b.rateTime(j));
    if (e <= s + 42.0 * QL_EPSILON)
        return 0.0;
    std::set<Time> knots;
    for (const Time k : a.integrationKnots(s, e))
        knots.insert(k);
    for (const Time k : b.integrationKnots(s, e))
        knots.insert(k);
    std::vector<Time> pieces(knots.begin(), knots.end());
    pieces.insert(pieces.begin(), s);
    pieces.push_back(e);
    static const GaussLegendreIntegration quad(8);
    Real sum = 0.0;
    for (Size p = 0; p + 1 < pieces.size(); ++p) {
        const Time lo = pieces[p], hi = pieces[p + 1];
        if (hi <= lo + tEps)
            continue;
        const Real half = 0.5 * (hi - lo), c = 0.5 * (lo + hi);
        sum += half * quad([&](const Real x) {
            const Time u = c + half * x;
            const Real fa = a.volLevel(i, u) * (decayI ? a.decay(i, u) : 1.0);
            const Real fb = b.volLevel(j, u) * (decayJ ? b.decay(j, u) : 1.0);
            return fa * fb;
        });
    }
    return sum;
}

Real FmmJointModel::volIntegral(const FmmParametrization& a, const Size i, const bool decayI, const Time s, Time e) {
    if (decayI)
        e = std::min(e, a.rateTime(i));
    if (e <= s + 42.0 * QL_EPSILON)
        return 0.0;
    std::vector<Time> pieces = a.integrationKnots(s, e);
    pieces.insert(pieces.begin(), s);
    pieces.push_back(e);
    static const GaussLegendreIntegration quad(8);
    Real sum = 0.0;
    for (Size p = 0; p + 1 < pieces.size(); ++p) {
        const Time lo = pieces[p], hi = pieces[p + 1];
        if (hi <= lo + tEps)
            continue;
        const Real half = 0.5 * (hi - lo), c = 0.5 * (lo + hi);
        sum += half * quad([&](const Real x) {
            const Time u = c + half * x;
            return a.volLevel(i, u) * (decayI ? a.decay(i, u) : 1.0);
        });
    }
    return sum;
}

FmmJointModel::StepData FmmJointModel::makeStep(const Time s, const Time e) const {
    QL_REQUIRE(e > s + tEps, "FmmJointModel::makeStep: need e > s");
    for (const Time bp : fxVolTimes_)
        QL_REQUIRE(bp <= s + tEps || bp >= e - tEps,
                   "FmmJointModel::makeStep: step [" << s << "," << e << "] crosses the FX vol breakpoint " << bp);
    StepData st;
    st.s = s;
    st.e = e;
    st.d = d_->makeStep(s, e);
    st.f = f_->makeStep(s, e);
    const Real sigmaX = fxVol(0.5 * (s + e));
    st.fxVariance = sigmaX * sigmaX * (e - s);
    const auto& pd = *d_->parametrization();
    const auto& pf = *f_->parametrization();
    const Size nd = Md_ + 1, nf = Mf_ + 1, n = nd + nf + 1;
    // the models' own shock covariances (rates and the live period's front stub), reconstructed
    // from their pseudo square roots so that the joint factorisation reproduces them exactly
    const Matrix covD = st.d.shockSqrt * transpose(st.d.shockSqrt);
    const Matrix covF = st.f.shockSqrt * transpose(st.f.shockSqrt);
    Matrix full(n, n, 0.0);
    for (Size a = 0; a < nd; ++a)
        for (Size b = 0; b < nd; ++b)
            full[a][b] = covD[a][b];
    for (Size a = 0; a < nf; ++a)
        for (Size b = 0; b < nf; ++b)
            full[nd + a][nd + b] = covF[a][b];
    // cross-currency block: level factors only. A rate's shock is lambda g times its loading on
    // the level factor L_{i1}; the front-stub shock (row M) is the live rate's level without decay.
    auto loadD = [&](const Size i) { return pd.factorLoadings()[i - 1][0]; };
    auto loadF = [&](const Size j) { return pf.factorLoadings()[j - 1][0]; };
    auto rowD = [&](const Size a, Real& load, Size& idx, bool& decay) {
        if (a < Md_) {
            idx = a + 1;
            decay = true;
        } else {
            idx = st.d.live;
            decay = false;
        }
        load = idx > 0 ? loadD(idx) : 0.0;
    };
    auto rowF = [&](const Size a, Real& load, Size& idx, bool& decay) {
        if (a < Mf_) {
            idx = a + 1;
            decay = true;
        } else {
            idx = st.f.live;
            decay = false;
        }
        load = idx > 0 ? loadF(idx) : 0.0;
    };
    for (Size a = 0; a < nd; ++a) {
        Real la;
        Size ia;
        bool da;
        rowD(a, la, ia, da);
        if (ia == 0)
            continue;
        for (Size b = 0; b < nf; ++b) {
            Real lb;
            Size ib;
            bool db;
            rowF(b, lb, ib, db);
            if (ib == 0)
                continue;
            const Real c = corr_.ratesDF * la * lb * crossIntegral(pd, ia, da, pf, ib, db, s, e);
            full[a][nd + b] = full[nd + b][a] = c;
        }
    }
    // FX row: sigma_X times the rates' vol integrals times the correlations with the level factors
    const Size x = nd + nf;
    full[x][x] = st.fxVariance;
    for (Size a = 0; a < nd; ++a) {
        Real la;
        Size ia;
        bool da;
        rowD(a, la, ia, da);
        if (ia == 0)
            continue;
        const Real c = corr_.fxD * la * sigmaX * volIntegral(pd, ia, da, s, e);
        full[a][x] = full[x][a] = c;
    }
    st.quantoDrift = Array(Mf_, 0.0);
    for (Size b = 0; b < nf; ++b) {
        Real lb;
        Size ib;
        bool db;
        rowF(b, lb, ib, db);
        if (ib == 0)
            continue;
        const Real integral = volIntegral(pf, ib, db, s, e);
        const Real c = corr_.fxF * lb * sigmaX * integral;
        full[nd + b][x] = full[x][nd + b] = c;
        // quanto drift of the foreign rates under the domestic spot measure (rates only, not the stub)
        if (b < Mf_)
            st.quantoDrift[b] = -c;
    }
    st.shockSqrt = pseudoSqrt(full, SalvagingAlgorithm::Spectral);
    return st;
}

void FmmJointModel::evolve(State& state, const StepData& step, const Array& z) const {
    QL_REQUIRE(z.size() == dimension(), "FmmJointModel::evolve: need " << dimension() << " draws, got " << z.size());
    QL_REQUIRE(std::fabs(state.t - step.s) < 1e-8, "FmmJointModel::evolve: state time " << state.t
                                                                                        << " != step start " << step.s);
    const Size nd = Md_ + 1, nf = Mf_ + 1, n = nd + nf + 1;
    Array v(n, 0.0);
    for (Size a = 0; a < n; ++a) {
        Real acc = 0.0;
        for (Size b = 0; b < n; ++b)
            acc += step.shockSqrt[a][b] * z[b];
        v[a] = acc;
    }
    Array vd(nd), vf(nf);
    for (Size a = 0; a < nd; ++a)
        vd[a] = v[a];
    for (Size a = 0; a < nf; ++a)
        vf[a] = v[nd + a];
    const Real bdBefore = d_->bankAccount(state.d), bfBefore = f_->bankAccount(state.f);
    d_->evolveWithCorrelatedShocks(state.d, step.d, vd);
    f_->evolveWithCorrelatedShocks(state.f, step.f, vf, step.quantoDrift);
    const Real bdAfter = d_->bankAccount(state.d), bfAfter = f_->bankAccount(state.f);
    // X B_f / B_d is a martingale under the domestic spot measure (exact discrete recursion)
    state.fx *= (bdAfter / bdBefore) / (bfAfter / bfBefore) * std::exp(v[n - 1] - 0.5 * step.fxVariance);
    state.t = step.e;
}

FmmJointPathGenerator::FmmJointPathGenerator(const QuantLib::ext::shared_ptr<FmmJointModel>& model,
                                             const std::vector<Time>& userTimes, const SequenceType sequenceType,
                                             const BigNatural seed, const SobolBrownianGenerator::Ordering ordering,
                                             const SobolRsg::DirectionIntegers directionIntegers)
    : model_(model), userTimes_(userTimes), sequenceType_(sequenceType), seed_(seed), ordering_(ordering),
      directionIntegers_(directionIntegers) {
    QL_REQUIRE(model_, "FmmJointPathGenerator: model is null");
    QL_REQUIRE(!userTimes_.empty(), "FmmJointPathGenerator: no simulation times");
    for (Size i = 0; i < userTimes_.size(); ++i)
        QL_REQUIRE(userTimes_[i] > (i == 0 ? 0.0 : userTimes_[i - 1]) + tEps,
                   "FmmJointPathGenerator: user times must be positive and strictly increasing");
    const auto& p = model_->domestic()->parametrization();
    std::set<Time> grid(userTimes_.begin(), userTimes_.end());
    const Time horizon = userTimes_.back();
    for (Size j = 1; j <= p->numberOfRates(); ++j)
        if (p->rateTime(j) < horizon - tEps)
            grid.insert(p->rateTime(j));
    for (const Time bp : model_->fxVolTimes())
        if (bp < horizon - tEps && bp > tEps)
            grid.insert(bp);
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

void FmmJointPathGenerator::reset() {
    variates_ = makeMultiPathVariateGenerator(sequenceType_, model_->dimension(), steps_.size(), seed_, ordering_,
                                              directionIntegers_);
}

FmmJointPathGenerator::Path FmmJointPathGenerator::next() const {
    const auto sample = variates_->next();
    QL_REQUIRE(sample.value.size() == steps_.size(), "FmmJointPathGenerator: unexpected variate dimensioning");
    Path path;
    FmmJointModel::State state = model_->initialState();
    Size nextUser = 0;
    for (Size i = 0; i < steps_.size(); ++i) {
        model_->evolve(state, steps_[i], sample.value[i]);
        while (nextUser < userIndex_.size() && userIndex_[nextUser] == i) {
            path.states.push_back(state);
            ++nextUser;
        }
    }
    QL_REQUIRE(path.states.size() == userTimes_.size(), "FmmJointPathGenerator: user time mapping failed");
    return path;
}

} // namespace QuantExt
