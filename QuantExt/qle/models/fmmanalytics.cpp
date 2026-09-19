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

#include <ql/pricingengines/blackformula.hpp>

namespace QuantExt {

namespace {
Real discount(const FmmParametrization& p, const Size m) {
    return p.termStructure()->discount(p.rateTime(m));
}
Real initialRate(const FmmParametrization& p, const Size j) {
    return (discount(p, j - 1) / discount(p, j) - 1.0) / p.tau(j);
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
    return (discount(p, swap.a) - discount(p, swap.b)) / fmmAnnuity(p, swap);
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
                                          const Option::Type type) {
    swap.validate(p.numberOfRates());
    FmmSwaptionApproxResult res;
    const Real Pa = discount(p, swap.a), Pb = discount(p, swap.b);
    const Real A = fmmAnnuity(p, swap);
    const Real S0 = (Pa - Pb) / A;
    res.forward = S0;
    res.annuity = A;

    // exact gradient q_i = dS/dR_i at t = 0. With c_i = tau_i/(1 + tau_i R_i(0)) and
    // A_i = sum of fixed-leg PVs paying at or after T_i:
    //   q_i = c_i [ (Pb 1_{i<=b} - Pa 1_{i<=a}) / A + S0 A_i / A ]
    // which vanishes identically for i <= a (front discounting cancels in S).
    const Size b = swap.b;
    std::vector<Real> q(b + 1, 0.0);
    for (Size i = swap.a + 1; i <= b; ++i) {
        const Real ci = p.tau(i) / (1.0 + p.tau(i) * initialRate(p, i));
        Real Ai = 0.0;
        for (Size c = 0; c < swap.fixedPayIndices.size(); ++c)
            if (swap.fixedPayIndices[c] >= i)
                Ai += swap.fixedAccruals[c] * discount(p, swap.fixedPayIndices[c]);
        q[i] = ci * (Pb / A + S0 * Ai / A);
    }

    // total normal variance with local vols frozen at time-0 rates; time integrals exact and
    // decay-aware up to the expiry T_a (the rho^eff factor sits inside integratedCovariance)
    const Time Te = p.rateTime(swap.a);
    QL_REQUIRE(Te > 0.0, "fmmSwaptionApprox: zero expiry");
    Real variance = 0.0;
    for (Size i = swap.a + 1; i <= b; ++i) {
        const Real phiI = p.phi(i, initialRate(p, i));
        for (Size k = swap.a + 1; k <= b; ++k) {
            const Real phiK = p.phi(k, initialRate(p, k));
            variance += q[i] * q[k] * phiI * phiK * p.integratedCovariance(i, k, 0.0, Te);
        }
    }
    res.variance = variance;
    res.normalVol = std::sqrt(variance / Te);
    res.price = A * bachelierBlackFormula(type, K, S0, std::sqrt(variance));
    return res;
}

} // namespace QuantExt
