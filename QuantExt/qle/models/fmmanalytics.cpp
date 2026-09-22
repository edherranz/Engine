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

#include <qle/models/fmmanalytics.hpp>
#include <qle/models/forwardmarketmodel.hpp>

#include <ql/math/statistics/incrementalstatistics.hpp>
#include <ql/pricingengines/blackformula.hpp>

namespace QuantExt {

namespace {
Real discount(const FmmParametrization& p, const Size m) {
    return p.termStructure()->discount(p.rateTime(m));
}
Real initialRate(const FmmParametrization& p, const Size j) {
    return (discount(p, j - 1) / discount(p, j) - 1.0) / p.tau(j);
}
// time-0 value of the deterministic basis amounts of the floating leg (per unit notional)
Real basisValue(const FmmParametrization& p, const FmmSwapSpec& swap) {
    Real v = 0.0;
    for (Size c = 0; c < swap.basisPayIndices.size(); ++c)
        v += swap.basisAmounts[c] * discount(p, swap.basisPayIndices[c]);
    return v;
}

// time-0 geometry shared by the approximations: annuity, forward, exact gradient q_i = dS/dR_i,
// local-vol factors phi_i at time-0 rates and the integrated covariances IC_ik over [0, T_a]
struct SwaptionGeometry {
    Real A = 0.0, S0 = 0.0;
    Time Te = 0.0;
    Size a = 0, b = 0;
    std::vector<Real> q, phi;          // indexed by rate index, entries outside (a, b] zero
    std::vector<std::vector<Real>> IC; // IC[i][k]
};

SwaptionGeometry geometry(const FmmParametrization& p, const FmmSwapSpec& swap) {
    SwaptionGeometry g;
    g.a = swap.a;
    g.b = swap.b;
    const Real Pa = discount(p, swap.a), Pb = discount(p, swap.b);
    g.A = fmmAnnuity(p, swap);
    g.S0 = (Pa - Pb + basisValue(p, swap)) / g.A;
    // exact gradient q_i = dS/dR_i at t = 0. With c_i = tau_i/(1 + tau_i R_i(0)),
    // A_i = sum of fixed-leg PVs paying at or after T_i and B_i = sum of basis PVs paying at or
    // after T_i (dP_j/dR_i = -c_i P_j for i <= j):
    //   q_i = c_i [ (Pb 1_{i<=b} - Pa 1_{i<=a} - B_i) / A + S0 A_i / A ]
    // which vanishes identically for i <= a (front discounting cancels in S).
    g.q.assign(g.b + 1, 0.0);
    g.phi.assign(g.b + 1, 0.0);
    for (Size i = swap.a + 1; i <= g.b; ++i) {
        const Real ci = p.tau(i) / (1.0 + p.tau(i) * initialRate(p, i));
        Real Ai = 0.0, Bi = 0.0;
        for (Size c = 0; c < swap.fixedPayIndices.size(); ++c)
            if (swap.fixedPayIndices[c] >= i)
                Ai += swap.fixedAccruals[c] * discount(p, swap.fixedPayIndices[c]);
        for (Size c = 0; c < swap.basisPayIndices.size(); ++c)
            if (swap.basisPayIndices[c] >= i)
                Bi += swap.basisAmounts[c] * discount(p, swap.basisPayIndices[c]);
        g.q[i] = ci * ((Pb - Bi) / g.A + g.S0 * Ai / g.A);
        g.phi[i] = p.phi(i, initialRate(p, i));
    }
    // time integrals exact and decay-aware up to the expiry T_a (the rho^eff factor sits inside
    // integratedCovariance)
    g.Te = p.rateTime(swap.a);
    QL_REQUIRE(g.Te > 0.0, "fmmSwaptionApprox: zero expiry");
    g.IC.assign(g.b + 1, std::vector<Real>(g.b + 1, 0.0));
    for (Size i = swap.a + 1; i <= g.b; ++i)
        for (Size k = swap.a + 1; k <= g.b; ++k)
            g.IC[i][k] = p.integratedCovariance(i, k, 0.0, g.Te);
    return g;
}
} // namespace

void FmmSwapSpec::validate(const Size M) const {
    QL_REQUIRE(a < b && b <= M, "FmmSwapSpec: need a < b <= M, got a=" << a << " b=" << b << " M=" << M);
    QL_REQUIRE(!fixedPayIndices.empty() && fixedPayIndices.size() == fixedAccruals.size(),
               "FmmSwapSpec: inconsistent fixed leg");
    for (Size c = 0; c < fixedPayIndices.size(); ++c) {
        QL_REQUIRE(fixedPayIndices[c] > a && fixedPayIndices[c] <= b,
                   "FmmSwapSpec: fixed pay index " << fixedPayIndices[c] << " outside (a, b]");
        QL_REQUIRE(c == 0 || fixedPayIndices[c] > fixedPayIndices[c - 1],
                   "FmmSwapSpec: fixed pay indices must be strictly increasing");
        QL_REQUIRE(fixedAccruals[c] > 0.0, "FmmSwapSpec: non-positive fixed accrual");
    }
    QL_REQUIRE(basisPayIndices.size() == basisAmounts.size(), "FmmSwapSpec: inconsistent basis flows");
    for (Size c = 0; c < basisPayIndices.size(); ++c) {
        QL_REQUIRE(basisPayIndices[c] > a && basisPayIndices[c] <= b,
                   "FmmSwapSpec: basis pay index " << basisPayIndices[c] << " outside (a, b]");
        QL_REQUIRE(c == 0 || basisPayIndices[c] > basisPayIndices[c - 1],
                   "FmmSwapSpec: basis pay indices must be strictly increasing");
    }
}

Real fmmCapletStdDev(const FmmParametrization& p, const Size j, const bool backwardLooking) {
    const Time Te = backwardLooking ? p.rateTime(j) : p.rateTime(j - 1);
    if (Te <= 0.0)
        return 0.0;
    // int_0^Te lambda_j^2 g_j^2, exact (the diagonal effective correlation is 1 by construction)
    return std::sqrt(p.integratedCovariance(j, j, 0.0, Te) / p.correlation(j, j));
}

Real fmmCapletPrice(const FmmParametrization& p, const Size j, const Real K, const bool backwardLooking,
                    const Option::Type type) {
    const Real F0 = initialRate(p, j);
    const Real sd = fmmCapletStdDev(p, j, backwardLooking);
    const Real df = discount(p, j);
    const Real tau = p.tau(j);
    if (p.volType() == FmmParametrization::LocalVolType::Normal)
        return df * tau * bachelierBlackFormula(type, K, F0, sd);
    return df * tau * blackFormula(type, K, F0, sd, 1.0, p.shift(j));
}

Real fmmCapletNormalVol(const FmmParametrization& p, const Size j, const Real K, const bool backwardLooking) {
    const Time Te = backwardLooking ? p.rateTime(j) : p.rateTime(j - 1);
    QL_REQUIRE(Te > 0.0, "fmmCapletNormalVol: zero expiry for rate " << j);
    if (p.volType() == FmmParametrization::LocalVolType::Normal)
        return fmmCapletStdDev(p, j, backwardLooking) / std::sqrt(Te);
    const Real F0 = initialRate(p, j);
    const Real undiscounted = fmmCapletPrice(p, j, K, backwardLooking) / (discount(p, j) * p.tau(j));
    return bachelierBlackFormulaImpliedVol(Option::Call, K, F0, Te, undiscounted);
}

Real fmmForwardSwapRate(const FmmParametrization& p, const FmmSwapSpec& swap) {
    swap.validate(p.numberOfRates());
    return (discount(p, swap.a) - discount(p, swap.b) + basisValue(p, swap)) / fmmAnnuity(p, swap);
}

Real fmmAnnuity(const FmmParametrization& p, const FmmSwapSpec& swap) {
    swap.validate(p.numberOfRates());
    Real annuity = 0.0;
    for (Size c = 0; c < swap.fixedPayIndices.size(); ++c)
        annuity += swap.fixedAccruals[c] * discount(p, swap.fixedPayIndices[c]);
    QL_REQUIRE(annuity > 0.0, "fmmAnnuity: non-positive annuity");
    return annuity;
}

FmmSwaptionApproxResult fmmSwaptionApprox(const FmmParametrization& p, const FmmSwapSpec& swap, const Real K,
                                          const Option::Type type, const FmmSwaptionApproxMethod method) {
    swap.validate(p.numberOfRates());
    const SwaptionGeometry g = geometry(p, swap);
    FmmSwaptionApproxResult res;
    res.forward = g.S0;
    res.annuity = g.A;
    res.method = method;

    // total normal variance with local vols frozen at time-0 rates
    Real V0 = 0.0;
    for (Size i = g.a + 1; i <= g.b; ++i)
        for (Size k = g.a + 1; k <= g.b; ++k)
            V0 += g.q[i] * g.q[k] * g.phi[i] * g.phi[k] * g.IC[i][k];
    QL_REQUIRE(V0 > 0.0, "fmmSwaptionApprox: non-positive swap-rate variance");

    auto normalResult = [&]() {
        res.variance = V0;
        res.normalVol = std::sqrt(V0 / g.Te);
        res.price = g.A * bachelierBlackFormula(type, K, g.S0, std::sqrt(V0));
        res.effectiveShift = Null<Real>();
        res.shiftedStdDev = 0.0;
        return res;
    };
    if (method == FmmSwaptionApproxMethod::StrikeFlatNormal || p.volType() == FmmParametrization::LocalVolType::Normal)
        return normalResult();

    // effective shift: project the local-vol factors onto the swap rate, beta_i = Cov(dR_i, dS) /
    // Var(dS) with the same integrated covariances; dphi_i/dR_i = 1 under displaced diffusion, so
    // the local variance slope at S(0) is 2 D with D = sum_ik q_i q_k beta_i phi_k IC_ik, and a
    // displaced-lognormal form V(S) = c^2 (S + delta_S)^2 matching V(0) and the slope has
    // 1/(S(0) + delta_S) = D / V(0)
    std::vector<Real> cov(g.b + 1, 0.0);
    for (Size i = g.a + 1; i <= g.b; ++i)
        for (Size k = g.a + 1; k <= g.b; ++k)
            cov[i] += g.q[k] * g.phi[i] * g.phi[k] * g.IC[i][k];
    Real D = 0.0;
    for (Size i = g.a + 1; i <= g.b; ++i) {
        const Real betaI = cov[i] / V0;
        for (Size k = g.a + 1; k <= g.b; ++k)
            D += g.q[i] * g.q[k] * betaI * g.phi[k] * g.IC[i][k];
    }
    const Real kappa = D / V0;
    if (kappa < 1e-8)
        return normalResult(); // no projected skew: normal limit
    const Real deltaS = 1.0 / kappa - g.S0;
    QL_REQUIRE(g.S0 + deltaS > 0.0 && K + deltaS > 0.0,
               "fmmSwaptionApprox: strike " << K << " outside the effective-shift domain (effective shift "
                                            << deltaS << ", forward " << g.S0 << ")");
    const Real sd = std::sqrt(V0) / (g.S0 + deltaS);
    res.effectiveShift = deltaS;
    res.shiftedStdDev = sd;
    res.price = g.A * blackFormula(type, K, g.S0, sd, 1.0, deltaS);
    // implied normal vol at the strike; the inversion is ill-conditioned when the time value is
    // at rounding level (tiny vol probed by a bootstrap solver, or a deep in-the-money strike),
    // where the exact small-vol limit sigma_N = sd (S0 + delta_S) / sqrt(T) applies
    const Real intrinsic = std::max((type == Option::Call ? 1.0 : -1.0) * (g.S0 - K), 0.0);
    const Real timeValue = res.price / g.A - intrinsic;
    if (timeValue > 1e-12 * std::max(1.0, intrinsic))
        res.normalVol = bachelierBlackFormulaImpliedVol(type, K, g.S0, g.Te, res.price / g.A);
    else
        res.normalVol = sd * (g.S0 + deltaS) / std::sqrt(g.Te);
    res.variance = res.normalVol * res.normalVol * g.Te;
    return res;
}

std::vector<FmmSwaptionMcResult> fmmSwaptionMc(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                                               const std::vector<FmmSwapSpec>& swaps, const std::vector<Real>& strikes,
                                               const Option::Type type, const Size pathsPerRep, const Size reps,
                                               const BigNatural seed) {
    QL_REQUIRE(model, "fmmSwaptionMc: model is null");
    QL_REQUIRE(!swaps.empty() && swaps.size() == strikes.size(), "fmmSwaptionMc: swaps/strikes mismatch");
    QL_REQUIRE(reps >= 2 && pathsPerRep > 0, "fmmSwaptionMc: need at least 2 replications");
    const auto& p = *model->parametrization();
    const Size a = swaps.front().a;
    for (const auto& s : swaps) {
        s.validate(p.numberOfRates());
        QL_REQUIRE(s.a == a, "fmmSwaptionMc: all swaps must share the expiry index");
    }
    const Time Te = p.rateTime(a);
    const Real w = type == Option::Call ? 1.0 : -1.0;
    std::vector<IncrementalStatistics> repMeans(swaps.size());
    for (Size r = 0; r < reps; ++r) {
        FmmPathGenerator gen(model, std::vector<Time>{Te}, Burley2020Sobol, seed + 13 * static_cast<BigNatural>(r));
        std::vector<IncrementalStatistics> acc(swaps.size());
        for (Size n = 0; n < pathsPerRep; ++n) {
            const auto path = gen.next();
            const auto& st = path.states[0];
            const Real bank = model->bankAccount(st);
            for (Size sw = 0; sw < swaps.size(); ++sw) {
                Real A = 0.0;
                for (Size c = 0; c < swaps[sw].fixedPayIndices.size(); ++c)
                    A += swaps[sw].fixedAccruals[c] * model->discountBond(st, p.rateTime(swaps[sw].fixedPayIndices[c]));
                Real basis = 0.0;
                for (Size c = 0; c < swaps[sw].basisPayIndices.size(); ++c)
                    basis += swaps[sw].basisAmounts[c] *
                             model->discountBond(st, p.rateTime(swaps[sw].basisPayIndices[c]));
                const Real S = (1.0 - model->discountBond(st, p.rateTime(swaps[sw].b)) + basis) / A;
                acc[sw].add(A * std::max(w * (S - strikes[sw]), 0.0) / bank);
            }
        }
        for (Size sw = 0; sw < swaps.size(); ++sw)
            repMeans[sw].add(acc[sw].mean());
    }
    std::vector<FmmSwaptionMcResult> out(swaps.size());
    for (Size sw = 0; sw < swaps.size(); ++sw) {
        const Real S0 = fmmForwardSwapRate(p, swaps[sw]);
        const Real A0 = fmmAnnuity(p, swaps[sw]);
        out[sw].price = repMeans[sw].mean();
        out[sw].priceSe = repMeans[sw].errorEstimate();
        out[sw].normalVol = bachelierBlackFormulaImpliedVol(type, strikes[sw], S0, Te, out[sw].price / A0);
        const Real volLo =
            bachelierBlackFormulaImpliedVol(type, strikes[sw], S0, Te, (out[sw].price - out[sw].priceSe) / A0);
        out[sw].normalVolSe = std::fabs(out[sw].normalVol - volLo);
    }
    return out;
}

} // namespace QuantExt
