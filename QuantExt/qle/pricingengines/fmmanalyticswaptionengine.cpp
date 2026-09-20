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

#include <qle/pricingengines/fmmanalyticswaptionengine.hpp>

#include <cmath>

namespace QuantExt {

FmmAnalyticSwaptionEngine::FmmAnalyticSwaptionEngine(const QuantLib::ext::shared_ptr<FmmParametrization>& parametrization,
                                                     const QuantLib::ext::shared_ptr<FmmGrid>& grid,
                                                     const FmmSwaptionApproxMethod method,
                                                     const Natural gridToleranceDays)
    : parametrization_(parametrization), grid_(grid), method_(method), gridToleranceDays_(gridToleranceDays) {
    QL_REQUIRE(parametrization_, "FmmAnalyticSwaptionEngine: parametrization is null");
    QL_REQUIRE(grid_, "FmmAnalyticSwaptionEngine: grid is null");
    QL_REQUIRE(grid_->numberOfRates() == parametrization_->numberOfRates(),
               "FmmAnalyticSwaptionEngine: grid (" << grid_->numberOfRates() << " periods) and parametrization ("
                                                   << parametrization_->numberOfRates() << ") disagree");
    registerWith(parametrization_->termStructure());
}

void FmmAnalyticSwaptionEngine::calculate() const {
    QL_REQUIRE(arguments_.settlementType == Settlement::Physical,
               "FmmAnalyticSwaptionEngine: physical settlement only (cash settlement needs an annuity convention)");
    Real strike = 0.0, nominal = 1.0;
    Option::Type type = Option::Call;
    const FmmSwapSpec spec = fmmSwapSpecFromSwaption(*grid_, arguments_, gridToleranceDays_, strike, type, nominal);
    // the swap's spread (if any) is not part of the FMM swap spec: treat a non-zero spread as an
    // adjustment of the fixed rate is NOT exact; reject explicitly
    QL_REQUIRE(arguments_.swap->spread() == 0.0 || std::fabs(arguments_.swap->spread()) < 1e-12,
               "FmmAnalyticSwaptionEngine: floating spreads are not supported in this release");
    const auto r = fmmSwaptionApprox(*parametrization_, spec, strike, type, method_);
    results_.value = nominal * r.price;
    results_.additionalResults["fmmNormalVol"] = r.normalVol;
    results_.additionalResults["fmmForwardSwapRate"] = r.forward;
    results_.additionalResults["fmmAnnuity"] = r.annuity;
    results_.additionalResults["fmmStrike"] = strike;
    results_.additionalResults["fmmExpiryTime"] = parametrization_->rateTime(spec.a);
    results_.additionalResults["fmmApproximationMethod"] =
        std::string(method_ == FmmSwaptionApproxMethod::EffectiveShift ? "EffectiveShift" : "StrikeFlatNormal");
    if (r.effectiveShift != Null<Real>())
        results_.additionalResults["fmmEffectiveShift"] = r.effectiveShift;
    results_.additionalResults["fmmGridMaxDateMismatchDays"] = static_cast<Real>(grid_->maxMismatchDays());
}

} // namespace QuantExt
