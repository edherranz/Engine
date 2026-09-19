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

#include <qle/models/fmmparametrization.hpp>
#include <qle/models/pseudoparameter.hpp>

#include <ql/math/comparison.hpp>
#include <ql/math/integrals/gaussianquadratures.hpp>
#include <ql/math/matrixutilities/symmetricschurdecomposition.hpp>

#include <algorithm>
#include <set>

namespace QuantExt {

FmmParametrization::FmmParametrization(const Currency& currency, const Handle<YieldTermStructure>& termStructure,
                                       const Array& rateTimes, const Array& shifts, const Array& volTimes,
                                       const std::vector<Array>& volLevels, const LocalVolType volType,
                                       const Real rhoInf, const Real beta, const Size factors, const std::string& name)
    : Parametrization(currency, name), termStructure_(termStructure), rateTimes_(rateTimes), shifts_(shifts),
      volTimes_(volTimes), volLevels_(volLevels), volType_(volType), M_(rateTimes.size() - 1), F_(factors) {
    checkGrid();
    QL_REQUIRE(shifts_.size() == M_, "FmmParametrization: shifts size (" << shifts_.size() << ") must equal M (" << M_
                                                                         << ")");
    QL_REQUIRE(volLevels_.size() == M_,
               "FmmParametrization: volLevels size (" << volLevels_.size() << ") must equal M (" << M_ << ")");
    for (Size j = 0; j < M_; ++j) {
        QL_REQUIRE(volLevels_[j].size() == volTimes_.size() + 1, "FmmParametrization: volLevels["
                                                                     << j << "] size (" << volLevels_[j].size()
                                                                     << ") must be volTimes size + 1 ("
                                                                     << volTimes_.size() + 1 << ")");
    }
    for (Size k = 1; k < volTimes_.size(); ++k)
        QL_REQUIRE(volTimes_[k] > volTimes_[k - 1], "FmmParametrization: volTimes must be strictly increasing");
    QL_REQUIRE(F_ >= 1 && F_ <= M_, "FmmParametrization: factors (" << F_ << ") must be in [1, " << M_ << "]");
    if (volType_ == LocalVolType::DisplacedDiffusion) {
        // admissibility (FMM_SPEC.md section 2.2): 0 <= delta_j <= 1/tau_j keeps the shifted-
        // lognormal support inside the positive-bond region 1 + tau_j R_j > 0
        for (Size j = 0; j < M_; ++j) {
            QL_REQUIRE(shifts_[j] >= 0.0, "FmmParametrization: negative shift for rate " << j + 1);
            QL_REQUIRE(shifts_[j] <= 1.0 / tau(j + 1) + 42.0 * QL_EPSILON,
                       "FmmParametrization: shift " << shifts_[j] << " for rate " << j + 1
                                                    << " exceeds the admissibility bound 1/tau = "
                                                    << 1.0 / tau(j + 1));
        }
    }
    makeFactorLoadings(rhoInf, beta);
    // raw parameter storage: concatenated levels, identity transform (calibration transforms in A3)
    volParameter_ = QuantLib::ext::make_shared<PseudoParameter>(M_ * (volTimes_.size() + 1));
    Size p = 0;
    for (Size j = 0; j < M_; ++j)
        for (Size k = 0; k < volLevels_[j].size(); ++k)
            volParameter_->setParam(p++, volLevels_[j][k]);
}

FmmParametrization::FmmParametrization(const Currency& currency, const Handle<YieldTermStructure>& termStructure,
                                       const Array& rateTimes,
                                       const QuantLib::ext::shared_ptr<IrLgm1fParametrization>& lgm,
                                       const std::string& name)
    : Parametrization(currency, name), termStructure_(termStructure), rateTimes_(rateTimes),
      volType_(LocalVolType::DisplacedDiffusion), M_(rateTimes.size() - 1), F_(1), lgm_(lgm) {
    checkGrid();
    QL_REQUIRE(lgm_, "FmmParametrization: LGM parametrization is null");
    shifts_ = Array(M_);
    for (Size j = 1; j <= M_; ++j)
        shifts_[j - 1] = 1.0 / tau(j);
    loadings_ = Matrix(M_, 1, 1.0);
    effectiveCorrelation_ = Matrix(M_, M_, 1.0);
}

void FmmParametrization::checkGrid() const {
    QL_REQUIRE(rateTimes_.size() >= 2, "FmmParametrization: need at least two rate times");
    QL_REQUIRE(rateTimes_[0] >= 0.0, "FmmParametrization: T_0 must be >= 0");
    for (Size k = 1; k < rateTimes_.size(); ++k)
        QL_REQUIRE(rateTimes_[k] > rateTimes_[k - 1] + QL_EPSILON,
                   "FmmParametrization: rate times must be strictly increasing");
    QL_REQUIRE(!termStructure_.empty(), "FmmParametrization: term structure is empty");
}

Time FmmParametrization::rateTime(const Size j) const {
    QL_REQUIRE(j <= M_, "FmmParametrization: rate time index " << j << " out of range [0," << M_ << "]");
    return rateTimes_[j];
}

Real FmmParametrization::tau(const Size j) const {
    QL_REQUIRE(j >= 1 && j <= M_, "FmmParametrization: tau index " << j << " out of range [1," << M_ << "]");
    // model-time accruals (FMM_SPEC.md assumption A8); schedule-accurate day counts arrive with
    // the OREData integration
    return rateTimes_[j] - rateTimes_[j - 1];
}

Real FmmParametrization::shift(const Size j) const {
    QL_REQUIRE(j >= 1 && j <= M_, "FmmParametrization: shift index " << j << " out of range [1," << M_ << "]");
    return shifts_[j - 1];
}

Size FmmParametrization::eta(const Time t) const {
    // smallest j >= 1 with T_j >= t; M+1 if beyond the grid
    Size j = 1;
    while (j <= M_ && rateTimes_[j] < t - 42.0 * QL_EPSILON)
        ++j;
    return j;
}

Real FmmParametrization::correlation(const Size i, const Size j) const {
    QL_REQUIRE(i >= 1 && i <= M_ && j >= 1 && j <= M_, "FmmParametrization: correlation index out of range");
    return effectiveCorrelation_[i - 1][j - 1];
}

Real FmmParametrization::volLevel(const Size j, const Time t) const {
    QL_REQUIRE(j >= 1 && j <= M_, "FmmParametrization: volLevel index " << j << " out of range [1," << M_ << "]");
    if (lgm_ != nullptr)
        return lgm_->alpha(t) * (lgm_->H(rateTimes_[j]) - lgm_->H(rateTimes_[j - 1]));
    Size k = std::upper_bound(volTimes_.begin(), volTimes_.end(), t - 42.0 * QL_EPSILON) - volTimes_.begin();
    return volLevels_[j - 1][k];
}

Real FmmParametrization::decay(const Size j, const Time t) const {
    QL_REQUIRE(j >= 1 && j <= M_, "FmmParametrization: decay index " << j << " out of range [1," << M_ << "]");
    const Time ts = rateTimes_[j - 1], te = rateTimes_[j];
    if (t <= ts)
        return 1.0;
    if (t >= te)
        return 0.0;
    if (lgm_ != nullptr) {
        const Real num = lgm_->H(te) - lgm_->H(t);
        const Real den = lgm_->H(te) - lgm_->H(ts);
        QL_REQUIRE(den > 0.0 || QuantLib::close_enough(den, 0.0),
                   "FmmParametrization: H must be non-decreasing for the implied decay");
        return QuantLib::close_enough(den, 0.0) ? 0.0 : std::min(1.0, std::max(0.0, num / den));
    }
    return (te - t) / (te - ts);
}

Real FmmParametrization::phi(const Size j, const Real r) const {
    return volType_ == LocalVolType::DisplacedDiffusion ? r + shift(j) : 1.0;
}

std::vector<Time> FmmParametrization::integrationKnots(const Time s, const Time e) const {
    std::set<Time> knots;
    for (Size k = 0; k < rateTimes_.size(); ++k)
        if (rateTimes_[k] > s + 42.0 * QL_EPSILON && rateTimes_[k] < e - 42.0 * QL_EPSILON)
            knots.insert(rateTimes_[k]);
    if (lgm_ == nullptr) {
        for (Size k = 0; k < volTimes_.size(); ++k)
            if (volTimes_[k] > s + 42.0 * QL_EPSILON && volTimes_[k] < e - 42.0 * QL_EPSILON)
                knots.insert(volTimes_[k]);
    } else {
        for (Size i = 0; i < lgm_->numberOfParameters(); ++i) {
            const Array& pt = lgm_->parameterTimes(i);
            for (Size k = 0; k < pt.size(); ++k)
                if (pt[k] > s + 42.0 * QL_EPSILON && pt[k] < e - 42.0 * QL_EPSILON)
                    knots.insert(pt[k]);
        }
    }
    return std::vector<Time>(knots.begin(), knots.end());
}

Real FmmParametrization::integratedCovariance(const Size i, const Size j, const Time s, Time e) const {
    QL_REQUIRE(i >= 1 && i <= M_ && j >= 1 && j <= M_, "FmmParametrization: covariance index out of range");
    // decays vanish after the earlier of the two maturities
    e = std::min(e, std::min(rateTimes_[i], rateTimes_[j]));
    if (e <= s + 42.0 * QL_EPSILON)
        return 0.0;
    std::vector<Time> knots = integrationKnots(s, e);
    knots.insert(knots.begin(), s);
    knots.push_back(e);
    Real sum = 0.0;
    for (Size k = 0; k + 1 < knots.size(); ++k)
        sum += intervalCovariance(i, j, knots[k], knots[k + 1], true, true);
    return sum * correlation(i, j);
}

Real FmmParametrization::integratedCovarianceSingleDecay(const Size i, const Size j, const Time s, Time e) const {
    QL_REQUIRE(i >= 1 && i <= M_ && j >= 1 && j <= M_, "FmmParametrization: covariance index out of range");
    e = std::min(e, rateTimes_[j]); // only the j-leg decays away
    if (e <= s + 42.0 * QL_EPSILON)
        return 0.0;
    std::vector<Time> knots = integrationKnots(s, e);
    knots.insert(knots.begin(), s);
    knots.push_back(e);
    Real sum = 0.0;
    for (Size k = 0; k + 1 < knots.size(); ++k)
        sum += intervalCovariance(i, j, knots[k], knots[k + 1], false, true);
    return sum * correlation(i, j);
}

Real FmmParametrization::integratedLevelVariance(const Size k, const Time s, const Time e) const {
    QL_REQUIRE(k >= 1 && k <= M_, "FmmParametrization: level variance index out of range");
    if (e <= s + 42.0 * QL_EPSILON)
        return 0.0;
    std::vector<Time> knots = integrationKnots(s, e);
    knots.insert(knots.begin(), s);
    knots.push_back(e);
    Real sum = 0.0;
    for (Size p = 0; p + 1 < knots.size(); ++p)
        sum += intervalCovariance(k, k, knots[p], knots[p + 1], false, false);
    return sum;
}

Real FmmParametrization::integratedDecay(const Size k, const Time s, Time e) const {
    QL_REQUIRE(k >= 1 && k <= M_, "FmmParametrization: decay integral index out of range");
    e = std::min(e, rateTimes_[k]);
    if (e <= s + 42.0 * QL_EPSILON)
        return 0.0;
    std::vector<Time> knots = integrationKnots(s, e);
    knots.insert(knots.begin(), s);
    knots.push_back(e);
    Real sum = 0.0;
    for (Size p = 0; p + 1 < knots.size(); ++p) {
        const Time a = knots[p], b = knots[p + 1];
        if (lgm_ != nullptr) {
            // g piecewise smooth in H; GL 16 per piece (exact for piecewise-linear H)
            static const GaussLegendreIntegration quad(16);
            const Real half = 0.5 * (b - a), c = 0.5 * (a + b);
            sum += half * quad([this, k, c, half](const Real x) { return decay(k, c + half * x); });
        } else {
            // linear decay: trapezoid is exact on a smooth piece
            sum += 0.5 * (decay(k, a) + decay(k, b)) * (b - a);
        }
    }
    return sum;
}

Real FmmParametrization::intervalCovariance(const Size i, const Size j, const Time s, const Time e,
                                            const bool decayOnI, const bool decayOnJ) const {
    // one smooth piece: no vol breakpoints or grid dates strictly inside (s, e)
    const Time mid = 0.5 * (s + e);
    if ((decayOnI && decay(i, mid) < QL_EPSILON) || (decayOnJ && decay(j, mid) < QL_EPSILON))
        return 0.0;

    if (lgm_ != nullptr) {
        // levels alpha(u) DeltaH; decayed legs replace DeltaH by (H(T) - H(u v T_{prev})); alpha
        // constant on the piece; exact for piecewise-linear H, machine precision otherwise (GL 16)
        static const GaussLegendreIntegration quad(16);
        const Real a2 = lgm_->alpha(mid) * lgm_->alpha(mid);
        const Real Hi = lgm_->H(rateTimes_[i]), Hj = lgm_->H(rateTimes_[j]);
        const Real dHi = Hi - lgm_->H(rateTimes_[i - 1]), dHj = Hj - lgm_->H(rateTimes_[j - 1]);
        const Time si = rateTimes_[i - 1], sj = rateTimes_[j - 1];
        auto f = [this, a2, Hi, Hj, dHi, dHj, si, sj, decayOnI, decayOnJ](const Real u) {
            const Real fi = decayOnI ? Hi - lgm_->H(std::max(u, si)) : dHi;
            const Real fj = decayOnJ ? Hj - lgm_->H(std::max(u, sj)) : dHj;
            return a2 * fi * fj;
        };
        const Real half = 0.5 * (e - s), c = 0.5 * (s + e);
        return half * quad([&f, half, c](const Real x) { return f(c + half * x); });
    }

    // piecewise-constant levels, linear decay: integrand is lambda_i lambda_j (A_i + B_i u)(A_j + B_j u)
    auto linearDecay = [this, mid](const Size idx, const bool withDecay, Real& A, Real& B) {
        const Time ts = rateTimes_[idx - 1], te = rateTimes_[idx];
        if (!withDecay || mid <= ts) {
            A = 1.0;
            B = 0.0;
        } else {
            A = te / (te - ts);
            B = -1.0 / (te - ts);
        }
    };
    Real Ai, Bi, Aj, Bj;
    linearDecay(i, decayOnI, Ai, Bi);
    linearDecay(j, decayOnJ, Aj, Bj);
    const Real li = volLevel(i, mid), lj = volLevel(j, mid);
    const Real d1 = e - s;
    const Real d2 = 0.5 * (e * e - s * s);
    const Real d3 = (e * e * e - s * s * s) / 3.0;
    return li * lj * (Ai * Aj * d1 + (Ai * Bj + Aj * Bi) * d2 + Bi * Bj * d3);
}

void FmmParametrization::makeFactorLoadings(const Real rhoInf, const Real beta) {
    QL_REQUIRE(rhoInf >= -1.0 && rhoInf <= 1.0, "FmmParametrization: rhoInf must be in [-1, 1]");
    QL_REQUIRE(beta >= 0.0, "FmmParametrization: beta must be >= 0");
    Matrix rho(M_, M_);
    for (Size a = 0; a < M_; ++a)
        for (Size b = 0; b < M_; ++b)
            rho[a][b] = rhoInf + (1.0 - rhoInf) * std::exp(-beta * std::fabs(rateTimes_[a + 1] - rateTimes_[b + 1]));
    SymmetricSchurDecomposition schur(rho);
    loadings_ = Matrix(M_, F_);
    for (Size a = 0; a < M_; ++a) {
        Real norm2 = 0.0;
        for (Size k = 0; k < F_; ++k) {
            const Real lambdaK = std::max(schur.eigenvalues()[k], 0.0);
            loadings_[a][k] = schur.eigenvectors()[a][k] * std::sqrt(lambdaK);
            norm2 += loadings_[a][k] * loadings_[a][k];
        }
        if (norm2 < QL_EPSILON) {
            loadings_[a][0] = 1.0;
            norm2 = 1.0;
        }
        const Real norm = std::sqrt(norm2);
        for (Size k = 0; k < F_; ++k)
            loadings_[a][k] /= norm;
    }
    effectiveCorrelation_ = loadings_ * transpose(loadings_);
}

const Array& FmmParametrization::parameterTimes(const Size i) const {
    QL_REQUIRE(lgm_ == nullptr && i == 0, "FmmParametrization: parameterTimes index out of range");
    return volTimes_;
}

const QuantLib::ext::shared_ptr<Parameter> FmmParametrization::parameter(const Size i) const {
    QL_REQUIRE(lgm_ == nullptr && i == 0, "FmmParametrization: parameter index out of range");
    return volParameter_;
}

void FmmParametrization::update() const {
    if (lgm_ != nullptr)
        return;
    // pull raw parameter values back into the level cache (identity transform)
    Size p = 0;
    for (Size j = 0; j < M_; ++j)
        for (Size k = 0; k < volLevels_[j].size(); ++k)
            volLevels_[j][k] = volParameter_->params()[p++];
}

} // namespace QuantExt
