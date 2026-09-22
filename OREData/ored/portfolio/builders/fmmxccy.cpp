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

#include <ored/model/fmmbuilder.hpp>
#include <ored/portfolio/builders/fmm.hpp>
#include <ored/utilities/dependencies.hpp>
#include <ored/utilities/indexparser.hpp>
#include <ored/utilities/log.hpp>
#include <ored/utilities/parsers.hpp>

#include <qle/models/fmmxccyreduction.hpp>
#include <qle/pricingengines/fmmlsmxccyreducedengine.hpp>
#include <qle/termstructures/correlationtermstructure.hpp>

#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>

#include <cmath>
#include <set>

using namespace QuantLib;
using namespace QuantExt;

namespace ore {
namespace data {

namespace {
Currency currencyOfKey(const std::string& key) {
    QuantLib::ext::shared_ptr<IborIndex> idx;
    if (tryParseIborIndex(key, idx))
        return idx->currency();
    return parseCurrency(key);
}
std::string rateFactorLabel(const std::string& key, const Currency& ccy) {
    // the name ORE's CrossAssetModel builders use for the currency's rate factor in the
    // correlation curves: the leg's index, or the currency's OIS index for a fixed leg
    QuantLib::ext::shared_ptr<IborIndex> idx;
    if (tryParseIborIndex(key, idx))
        return key;
    return swapIndexDiscountCurve(ccy.code());
}
} // namespace

/*! Hagan's reduced pricing of a cancellable cross-currency swap (docs/A9_1_HAGAN_REDUCTION_PLAN.md,
    A9_XCCY_TEMPLATE_PLAN.md section 5): the modelled currency B is the fixed leg's, the floating
    currency A is telescoped away; the B-currency FMM is calibrated, per exercise date, to the
    coterminal swaption at the strike R_fix - s_k^eq (eq. 3.3a) with the market normal vol
    (XccyMethod Standard) or the effective vol of eq. 3.41 (Improved, the default). The model
    curve is B's pricing discount curve (USD collateral: EUR-in-USD), the calibration helpers stay
    on the in-currency curves as for ORE's CrossAssetModel. */
QuantLib::ext::shared_ptr<PricingEngine> fmmXccyReducedEngine(const EngineBuilder* builder, const std::string& id,
                                                              const std::vector<std::string>& keys,
                                                              const std::vector<Date>& dates,
                                                              const std::vector<Date>& maturities,
                                                              const std::vector<std::vector<Real>>& strikes) {
    QL_REQUIRE(keys.size() == 2 && strikes.size() == 2,
               "FmmLsmSwaptionEngineBuilder: cross-currency swaptions need two currency keys, got " << keys.size());
    auto hasStrike = [](const std::vector<Real>& v) {
        for (const Real x : v)
            if (x != Null<Real>())
                return true;
        return false;
    };
    const bool s0 = hasStrike(strikes[0]), s1 = hasStrike(strikes[1]);
    QL_REQUIRE(s0 != s1, "FmmLsmSwaptionEngineBuilder: a cross-currency swaption needs exactly one fixed leg "
                         "(calibration strikes in one currency, an overnight floating leg in the other), keys "
                             << keys[0] << " / " << keys[1]);
    const Size iB = s0 ? 0 : 1, iA = 1 - iB;
    const Currency ccyA = currencyOfKey(keys[iA]), ccyB = currencyOfKey(keys[iB]);
    QL_REQUIRE(ccyA != ccyB, "FmmLsmSwaptionEngineBuilder: the two keys are in the same currency " << ccyA.code());
    const std::string method = builder->modelParameter("XccyMethod", {}, false, "Improved");
    QL_REQUIRE(method == "Standard" || method == "Improved",
               "FMM ModelParameter XccyMethod must be Standard or Improved, got " << method);
    const std::string pricingConfig = builder->configuration(MarketContext::pricing);
    auto market = builder->market();
    DLOG("FMM cross-currency reduction for " << id << ": B = " << keys[iB] << " (" << ccyB.code() << "), A = " << keys[iA]
                                             << " (" << ccyA.code() << "), method " << method);

    FmmModelRequest req;
    req.id = id;
    req.key = keys[iB];
    req.expiries = dates;
    req.maturities = maturities;
    req.strikes = strikes[iB];
    std::set<Date> gridDates(dates.begin(), dates.end());
    gridDates.insert(maturities.begin(), maturities.end());
    req.gridDates.assign(gridDates.begin(), gridDates.end());
    auto fmmBuilder = fmmModelBuilder(builder, req);

    // curves and spot of the reduction (collateral-consistent pricing curves); the model curve is
    // B's pricing discount curve, observed for recalibration (the base builder observes only the
    // calibration curve)
    Handle<YieldTermStructure> discountA = market->discountCurve(ccyA.code(), pricingConfig);
    Handle<YieldTermStructure> discountB = market->discountCurve(ccyB.code(), pricingConfig);
    Handle<Quote> fxSpotBperA = market->fxRate(ccyA.code() + ccyB.code(), pricingConfig);
    fmmBuilder->discountCurve().linkTo(*discountB);
    fmmBuilder->addObservable(*discountA);
    fmmBuilder->addObservable(*discountB);
    fmmBuilder->addObservable(*fxSpotBperA);
    Handle<SwaptionVolatilityStructure> svtsB = market->swaptionVol(keys[iB], pricingConfig);
    QL_REQUIRE(svtsB->volatilityType() == Normal,
               "FmmLsmSwaptionEngineBuilder: the " << keys[iB] << " swaption surface must be quoted in normal vols");
    Handle<SwaptionVolatilityStructure> svtsA;
    Handle<BlackVolTermStructure> fxVol;
    Handle<QuantExt::CorrelationTermStructure> rhoFxB, rhoAB;
    const std::string labelA = rateFactorLabel(keys[iA], ccyA), labelB = rateFactorLabel(keys[iB], ccyB);
    // Hagan's F^{A/B} is quoted in B per A; ORE's FX-GENERIC-A-B index is A's price in B, the
    // market's inverse pairs are handled (sign flip) by the correlation lookup
    const std::string fxIndexBperA = "FX-GENERIC-" + ccyA.code() + "-" + ccyB.code();
    if (method == "Improved") {
        svtsA = market->swaptionVol(keys[iA], pricingConfig);
        QL_REQUIRE(svtsA->volatilityType() == Normal,
                   "FmmLsmSwaptionEngineBuilder: the " << keys[iA] << " swaption surface must be quoted in normal vols");
        fxVol = market->fxVol(ccyB.code() + ccyA.code(), pricingConfig);
        rhoFxB = market->correlationCurve(fxIndexBperA, labelB, pricingConfig);
        rhoAB = market->correlationCurve(labelA, labelB, pricingConfig);
        fmmBuilder->addObservable(*svtsA);
        fmmBuilder->addObservable(*fxVol);
        fmmBuilder->addObservable(*rhoFxB);
        fmmBuilder->addObservable(*rhoAB);
    }

    FmmXccyReductionInputs inputs;
    inputs.ccyA = ccyA;
    inputs.ccyB = ccyB;
    inputs.discountA = discountA;
    inputs.discountB = discountB;
    inputs.fxSpotBperA = fxSpotBperA;
    const Currency npvCcy = currencyOfKey(keys[0]);

    auto gridProvider = [fmmBuilder](const std::set<Date>& d) {
        fmmBuilder->setGridDates(d);
        return fmmBuilder->prepareGrid();
    };
    auto modelProvider = [=](const FmmXccyReductionDiagnostics& diag, const FmmGrid& grid) {
        FmmXccyCalibrationRecord rec;
        rec.method = method;
        const std::vector<Date> basketExpiries = fmmBuilder->basketExpiryDates();
        std::vector<std::pair<Real, Real>> overrides;
        for (Size h = 0; h < basketExpiries.size(); ++h) {
            const FmmXccyReductionRight* r = nullptr;
            for (const auto& x : diag.rights)
                if (x.noticeDate == basketExpiries[h]) {
                    r = &x;
                    break;
                }
            QL_REQUIRE(r, "FmmLsmSwaptionEngineBuilder: no exercise right for the basket expiry "
                              << basketExpiries[h] << " (the basket must be the coterminal swaptions on the notice dates)");
            const Real sigmaB = svtsB->volatility(r->expiry, r->swapLengthB, r->strike);
            Real vol = sigmaB, sigmaA = 0.0, sigmaFx = 0.0, rFxB = 0.0, rAB = 0.0, omitted = 0.0;
            if (method == "Improved") {
                sigmaA = svtsA->volatility(r->expiry, r->swapLengthA, r->forwardSwapRateA);
                // ATM-forward FX option in the surface's own quotation (A per B)
                const Real fwdAperB =
                    discountB->discount(r->expiry) / discountA->discount(r->expiry) / fxSpotBperA->value();
                sigmaFx = fxVol->blackVol(r->expiry, fwdAperB);
                rFxB = rhoFxB->correlation(r->expiry);
                rAB = rhoAB->correlation(r->expiry);
                const Real s = r->equivalentSpread;
                // Hagan eq. 3.41 and the omitted-effect scale of eq. 3.46
                vol = sigmaB + rFxB * sigmaFx * s + r->gammaB * sigmaB * s - rAB * r->gammaA * sigmaA * s;
                omitted = (sigmaFx * sigmaFx + std::pow(r->gammaB * sigmaB, 2) + std::pow(r->gammaA * sigmaA, 2)) *
                          std::fabs(s);
                QL_REQUIRE(vol > 0.0, "FmmLsmSwaptionEngineBuilder: the effective volatility of the right at "
                                          << r->noticeDate << " is not positive (" << vol << "); no clipping");
            }
            overrides.push_back({r->strike, vol});
            rec.strikes.push_back(r->strike);
            rec.marketVolsB.push_back(sigmaB);
            rec.targetVols.push_back(vol);
            rec.sigmaA.push_back(sigmaA);
            rec.sigmaFx.push_back(sigmaFx);
            rec.rhoFxB.push_back(rFxB);
            rec.rhoAB.push_back(rAB);
            rec.omittedEffectVol.push_back(omitted);
        }
        fmmBuilder->setTargetOverrides(overrides);
        rec.model = fmmBuilder->modelAsFmm();
        QL_REQUIRE(fmmBuilder->grid()->numberOfRates() == grid.numberOfRates(),
                   "FmmLsmSwaptionEngineBuilder: the calibrated model's grid differs from the reduction grid");
        rec.extraResults["fmmXccyRateFactorA"] = labelA;
        rec.extraResults["fmmXccyRateFactorB"] = labelB;
        rec.extraResults["fmmXccyFxIndexBperA"] = fxIndexBperA;
        return rec;
    };
    FmmLsmEngineConfig cfg = fmmLsmEngineConfig(builder, fmmBuilder, false);
    return QuantLib::ext::make_shared<FmmLsmXccyReducedEngine>(inputs, npvCcy, gridProvider, modelProvider, cfg);
}

} // namespace data
} // namespace ore
