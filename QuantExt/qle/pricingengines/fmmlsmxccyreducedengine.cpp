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

#include <qle/pricingengines/fmmlsmxccyreducedengine.hpp>

#include <ql/settings.hpp>

#include <memory>

namespace QuantExt {

FmmLsmXccyReducedEngine::FmmLsmXccyReducedEngine(const FmmXccyReductionInputs& inputs, const Currency& npvCurrency,
                                                 const GridProvider& gridProvider, const ModelProvider& modelProvider,
                                                 const FmmLsmEngineConfig& config)
    : inputs_(inputs), npvCurrency_(npvCurrency), gridProvider_(gridProvider), modelProvider_(modelProvider),
      config_(config) {
    QL_REQUIRE(gridProvider_ && modelProvider_, "FmmLsmXccyReducedEngine: grid and model providers are required");
    QL_REQUIRE(!inputs_.discountA.empty() && !inputs_.discountB.empty() && !inputs_.fxSpotBperA.empty(),
               "FmmLsmXccyReducedEngine: discount curves and the FX spot are required");
    QL_REQUIRE(npvCurrency_ == inputs_.ccyA || npvCurrency_ == inputs_.ccyB,
               "FmmLsmXccyReducedEngine: the NPV currency " << npvCurrency_.code() << " is neither "
                                                             << inputs_.ccyA.code() << " nor " << inputs_.ccyB.code());
    registerWith(inputs_.discountA);
    registerWith(inputs_.discountB);
    registerWith(inputs_.fxSpotBperA);
    for (const auto& o : config_.observables)
        registerWith(o);
}

void FmmLsmXccyReducedEngine::calculate() const {
    const Date today = inputs_.discountB->referenceDate();
    QL_REQUIRE(inputs_.discountA->referenceDate() == today,
               "FmmLsmXccyReducedEngine: discount curve reference dates differ (" << inputs_.discountA->referenceDate()
                                                                                 << " vs " << today << ")");
    QL_REQUIRE(arguments_.legs.size() == arguments_.payer.size() && arguments_.legs.size() == arguments_.currency.size(),
               "FmmLsmXccyReducedEngine: inconsistent legs / payer / currency vectors");

    // 1. the grid on the trade's dates, 2. the reduction, 3. the model calibrated to the reduction
    const std::set<Date> dates = fmmXccyContractualDates(arguments_, today);
    auto grid = gridProvider_(dates);
    QL_REQUIRE(grid, "FmmLsmXccyReducedEngine: no grid");
    QL_REQUIRE(grid->referenceDate() == today,
               "FmmLsmXccyReducedEngine: grid reference date " << grid->referenceDate() << " != " << today);
    FmmXccyReductionDiagnostics diag;
    std::vector<Date> usedSettle;
    FmmCallableInstrument inst =
        fmmMapXccyReduced(arguments_, *grid, config_.gridToleranceDays, today, inputs_, diag, usedSettle);
    const FmmXccyCalibrationRecord rec = modelProvider_(diag, *grid);
    QL_REQUIRE(rec.model, "FmmLsmXccyReducedEngine: no model");
    const auto& p = *rec.model->parametrization();
    QL_REQUIRE(p.currency() == inputs_.ccyB, "FmmLsmXccyReducedEngine: model currency "
                                                 << p.currency().code() << " != " << inputs_.ccyB.code());
    QL_REQUIRE(p.numberOfRates() == grid->numberOfRates(),
               "FmmLsmXccyReducedEngine: model and grid disagree on the number of periods");
    QL_REQUIRE(p.termStructure()->referenceDate() == today,
               "FmmLsmXccyReducedEngine: model curve reference date " << p.termStructure()->referenceDate()
                                                                      << " != " << today);
    const Real toNpv = npvCurrency_ == inputs_.ccyB ? 1.0 : 1.0 / inputs_.fxSpotBperA->value();

    // reduction record (model currency B unless stated)
    auto& out = results_.additionalResults;
    // the reduced underlying: the flows after today plus a telescoped coupon start on the valuation
    // date itself (index 0, the forward-start value of the A leg), which the curve value skips
    const Real underlyingB = fmmUnderlyingCurveValue(inst, p) + inst.fixedFlows.front();
    results_.underlyingNpv = underlyingB * toNpv;
    out["underlyingNpv"] = underlyingB * toNpv;
    out["fmmXccyMethod"] = rec.method;
    out["fmmXccyCurrencyA"] = inputs_.ccyA.code();
    out["fmmXccyCurrencyB"] = inputs_.ccyB.code();
    out["fmmXccyNpvCurrency"] = npvCurrency_.code();
    out["fmmXccyFxSpotBperA"] = diag.fxSpotBperA;
    out["fmmXccyUnderlyingValueB"] = underlyingB;
    out["fmmXccySpreadStreamValueB"] = diag.spreadStreamValueB;
    out["fmmXccyTelescopedCoupons"] = static_cast<Real>(diag.telescopedCoupons);
    out["fmmXccyConvertedFlows"] = static_cast<Real>(diag.convertedFlows);
    std::vector<Date> noticeDates, settleDates;
    std::vector<Real> eqSpreads, strikes, fixedRates, annA, annB, gamA, gamB, fwdA, fwdFx, fees;
    for (const auto& r : diag.rights) {
        noticeDates.push_back(r.noticeDate);
        settleDates.push_back(r.settleDate);
        eqSpreads.push_back(r.equivalentSpread);
        strikes.push_back(r.strike);
        fixedRates.push_back(r.fixedRate);
        annA.push_back(r.annuityA);
        annB.push_back(r.annuityB);
        gamA.push_back(r.gammaA);
        gamB.push_back(r.gammaB);
        fwdA.push_back(r.forwardSwapRateA);
        fwdFx.push_back(r.forwardFxAtSettle);
        fees.push_back(r.feeFlowB);
    }
    out["fmmXccyNoticeDates"] = noticeDates;
    out["fmmXccySettlementDates"] = settleDates;
    out["fmmXccyEquivalentSpreads"] = eqSpreads;
    out["fmmXccyStrikes"] = strikes;
    out["fmmXccyFixedRates"] = fixedRates;
    out["fmmXccyAnnuityA"] = annA;
    out["fmmXccyAnnuityB"] = annB;
    out["fmmXccyGammaA"] = gamA;
    out["fmmXccyGammaB"] = gamB;
    out["fmmXccyForwardSwapRateA"] = fwdA;
    out["fmmXccyForwardFxAtSettle"] = fwdFx;
    out["fmmXccyFeeFlowsB"] = fees;
    out["fmmXccyCalibrationStrikes"] = rec.strikes;
    out["fmmXccyMarketVolsB"] = rec.marketVolsB;
    out["fmmXccyTargetVols"] = rec.targetVols;
    out["fmmXccyVolA"] = rec.sigmaA;
    out["fmmXccyVolFx"] = rec.sigmaFx;
    out["fmmXccyRhoFxB"] = rec.rhoFxB;
    out["fmmXccyRhoAB"] = rec.rhoAB;
    out["fmmXccyOmittedEffectVol"] = rec.omittedEffectVol;
    for (const auto& kv : rec.extraResults)
        out[kv.first] = kv.second;

    if (inst.rights.empty()) {
        // no exercise object: the instrument is the reduced underlying; no right left: worthless
        results_.value = arguments_.exercise ? 0.0 : underlyingB * toNpv;
        out["fmmGridMaxDateMismatchDays"] = static_cast<Real>(grid->maxMismatchDays());
        return;
    }

    FmmLsmPricer pricer(rec.model, inst, config_.lsm);
    FmmLsmResult res;
    bool reusedPolicy = false;
    if (config_.policyMode == FmmPolicyMode::Frozen && frozenPolicy_) {
        res = pricer.valueWithPolicy(*frozenPolicy_, config_.lsm.valuationSeed);
        reusedPolicy = true;
    } else {
        res = pricer.calculate();
        if (config_.policyMode == FmmPolicyMode::Frozen)
            frozenPolicy_ = pricer.policy();
    }
    results_.value = res.lowerBound * toNpv;
    std::unique_ptr<FmmDualBoundResult> dual;
    if (config_.dualBound && !reusedPolicy)
        dual = std::make_unique<FmmDualBoundResult>(
            pricer.dualBound(config_.dualOuterPaths, config_.dualInnerPaths, config_.dualSeed));
    fmmWriteLsmResults(out, res, dual.get(), config_, *grid); // LSM statistics in the model currency B
    out["fmmXccyValueB"] = res.lowerBound;
    out["fmmXccyValueNpvCurrency"] = res.lowerBound * toNpv;
    out["fmmExerciseDates"] = arguments_.exercise->dates();
    out["fmmSettlementDates"] = usedSettle;
    out["fmmPolicyMode"] =
        std::string(config_.policyMode == FmmPolicyMode::Frozen
                        ? (reusedPolicy ? "Frozen (policy of the first valuation reused)" : "Frozen (policy trained)")
                        : "Retrain");
    if (config_.calibrationResults)
        for (const auto& kv : config_.calibrationResults())
            out[kv.first] = kv.second;
}

} // namespace QuantExt
