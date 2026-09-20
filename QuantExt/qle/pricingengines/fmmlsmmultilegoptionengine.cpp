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

#include <qle/cashflows/overnightindexedcoupon.hpp>
#include <qle/pricingengines/fmmlsmmultilegoptionengine.hpp>

#include <ql/cashflows/capflooredcoupon.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/cashflows/iborcoupon.hpp>
#include <ql/cashflows/overnightindexedcoupon.hpp>
#include <ql/cashflows/simplecashflow.hpp>

#include <cmath>
#include <memory>
#include <sstream>

namespace QuantExt {

FmmLegMapper::FmmLegMapper(const FmmGrid& grid, const Natural toleranceDays, const Date& today)
    : grid_(grid), toleranceDays_(toleranceDays), today_(today) {}

void FmmLegMapper::addLeg(const Leg& leg, const Real sign, FmmCallableInstrument& inst, const std::string& what) const {
    const Size M = grid_.numberOfRates();
    if (inst.fixedFlows.size() != M + 1) {
        inst.fixedFlows.assign(M + 1, 0.0);
        inst.floatWeights.assign(M + 1, 0.0);
    }
    for (const auto& cf : leg) {
        if (cf->date() <= today_)
            continue; // already paid (or paying today): not part of the forward value
        if (auto c = QuantLib::ext::dynamic_pointer_cast<CappedFlooredCoupon>(cf)) {
            QL_FAIL("FmmLegMapper: " << what << ": capped/floored coupons are not supported in this release");
        } else if (auto on = QuantLib::ext::dynamic_pointer_cast<QuantExt::OvernightIndexedCoupon>(cf)) {
            QL_REQUIRE(on->accrualStartDate() >= today_,
                       "FmmLegMapper: " << what << ": coupon accruing at the valuation date (" << on->accrualStartDate()
                                        << " - " << on->accrualEndDate() << ") is not supported (already-started underlying)");
            QL_REQUIRE(!on->includeSpread(), "FmmLegMapper: " << what << ": includeSpread compounding is not supported");
            QL_REQUIRE(on->lookback() == 0 * Days && on->rateCutoff() == 0,
                       "FmmLegMapper: " << what
                                        << ": lookback / rate cut-off conventions are not supported in this "
                                           "release (baseline scope: plain compounded RFR)");
            FmmCallableInstrument::CompoundedFloat f;
            f.startIdx = grid_.index(on->accrualStartDate(), toleranceDays_, what + " accrual start");
            f.endIdx = grid_.index(on->accrualEndDate(), toleranceDays_, what + " accrual end");
            f.payIdx = grid_.index(on->date(), toleranceDays_, what + " pay date");
            QL_REQUIRE(f.endIdx > f.startIdx && f.payIdx >= f.endIdx, "FmmLegMapper: " << what << ": inconsistent coupon dates");
            f.weight = sign * on->nominal() * on->gearing();
            f.spreadAmount = sign * on->nominal() * on->spread() * on->accrualPeriod();
            inst.compoundedFloats.push_back(f);
            inst.lastFlowIdx = std::max(inst.lastFlowIdx, f.payIdx);
        } else if (auto onq = QuantLib::ext::dynamic_pointer_cast<QuantLib::OvernightIndexedCoupon>(cf)) {
            QL_REQUIRE(onq->accrualStartDate() >= today_,
                       "FmmLegMapper: " << what << ": coupon accruing at the valuation date is not supported");
            QL_REQUIRE(onq->lockoutDays() == 0,
                       "FmmLegMapper: " << what << ": lockout conventions are not supported in this release");
            FmmCallableInstrument::CompoundedFloat f;
            f.startIdx = grid_.index(onq->accrualStartDate(), toleranceDays_, what + " accrual start");
            f.endIdx = grid_.index(onq->accrualEndDate(), toleranceDays_, what + " accrual end");
            f.payIdx = grid_.index(onq->date(), toleranceDays_, what + " pay date");
            QL_REQUIRE(f.endIdx > f.startIdx && f.payIdx >= f.endIdx, "FmmLegMapper: " << what << ": inconsistent coupon dates");
            f.weight = sign * onq->nominal() * onq->gearing();
            f.spreadAmount = sign * onq->nominal() * onq->spread() * onq->accrualPeriod();
            inst.compoundedFloats.push_back(f);
            inst.lastFlowIdx = std::max(inst.lastFlowIdx, f.payIdx);
        } else if (auto ib = QuantLib::ext::dynamic_pointer_cast<IborCoupon>(cf)) {
            QL_FAIL("FmmLegMapper: " << what
                                     << ": term-rate (ibor) coupons are outside the single-currency RFR baseline of "
                                        "this release");
        } else if (auto fx = QuantLib::ext::dynamic_pointer_cast<FixedRateCoupon>(cf)) {
            QL_REQUIRE(fx->accrualStartDate() >= today_ || fx->accrualEndDate() == fx->date(),
                       "FmmLegMapper: " << what << ": fixed coupon accruing at the valuation date");
            const Size idx = grid_.index(fx->date(), toleranceDays_, what + " fixed pay date");
            inst.fixedFlows[idx] += sign * fx->amount();
            inst.lastFlowIdx = std::max(inst.lastFlowIdx, idx);
        } else if (QuantLib::ext::dynamic_pointer_cast<Coupon>(cf)) {
            QL_FAIL("FmmLegMapper: " << what << ": unsupported coupon type");
        } else {
            // notional exchanges, redemptions, fees
            const Size idx = grid_.index(cf->date(), toleranceDays_, what + " cash flow date");
            inst.fixedFlows[idx] += sign * cf->amount();
            inst.lastFlowIdx = std::max(inst.lastFlowIdx, idx);
        }
    }
}

Real FmmLegMapper::underlyingCurveValue(const FmmCallableInstrument& inst, const FmmParametrization& p) {
    return fmmUnderlyingCurveValue(inst, p);
}

void fmmWriteLsmResults(std::map<std::string, QuantLib::ext::any>& out, const FmmLsmResult& res,
                        const FmmDualBoundResult* dual, const FmmLsmEngineConfig& config, const FmmGrid& grid) {
    out["fmmLsmLowerBound"] = res.lowerBound;
    out["fmmLsmLowerBoundStdError"] = res.lowerBoundSe;
    out["fmmLsmTrainingValue"] = res.trainingValue;
    out["fmmLsmTrainingPaths"] = static_cast<Real>(config.lsm.trainingPaths);
    out["fmmLsmValuationPaths"] = static_cast<Real>(config.lsm.valuationPaths);
    out["fmmLsmTrainingSeed"] = static_cast<Real>(config.lsm.trainingSeed);
    out["fmmLsmValuationSeed"] = static_cast<Real>(config.lsm.valuationSeed);
    out["fmmLsmBasisOrder"] = static_cast<Real>(config.lsm.basisOrder);
    out["fmmNoExerciseProbability"] = res.noExerciseProbability;
    out["fmmExpectedExerciseTime"] = res.expectedExerciseTime;
    out["fmmExerciseProbabilityType"] = std::string("unconditional");
    std::vector<Real> probs = res.exerciseProbability, probSe = res.exerciseProbabilitySe;
    out["fmmExerciseProbabilities"] = probs;
    out["fmmExerciseProbabilityStdErrors"] = probSe;
    if (dual) {
        out["fmmDualityGap"] = dual->gap;
        out["fmmDualityGapStdError"] = dual->gapSe;
        out["fmmLsmUpperBound"] = res.lowerBound + dual->gap;
        out["fmmLsmUpperBoundStdError"] = std::sqrt(res.lowerBoundSe * res.lowerBoundSe + dual->gapSe * dual->gapSe);
        out["fmmDualOuterPaths"] = static_cast<Real>(dual->outerPaths);
        out["fmmDualInnerPaths"] = static_cast<Real>(dual->innerPaths);
    }
    out["fmmGridPeriods"] = static_cast<Real>(grid.numberOfRates());
    out["fmmGridMaxDateMismatchDays"] = static_cast<Real>(grid.maxMismatchDays());
    out["fmmRuntimeSeconds"] = res.runtimeSeconds + (dual ? dual->runtimeSeconds : 0.0);
}

FmmLsmMultiLegOptionEngine::FmmLsmMultiLegOptionEngine(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                                                       const QuantLib::ext::shared_ptr<FmmGrid>& grid,
                                                       const FmmLsmEngineConfig& config)
    : model_(model), grid_(grid), config_(config) {
    QL_REQUIRE(model_ && grid_, "FmmLsmMultiLegOptionEngine: model or grid is null");
    QL_REQUIRE(grid_->numberOfRates() == model_->parametrization()->numberOfRates(),
               "FmmLsmMultiLegOptionEngine: grid and model disagree on the number of periods");
    registerWith(model_->parametrization()->termStructure());
}

void FmmLsmMultiLegOptionEngine::calculate() const {
    const auto& p = *model_->parametrization();
    const Date today = grid_->referenceDate();
    QL_REQUIRE(p.termStructure()->referenceDate() == today,
               "FmmLsmMultiLegOptionEngine: model curve reference date " << p.termStructure()->referenceDate()
                                                                         << " != grid reference date " << today);
    QL_REQUIRE(arguments_.legs.size() == arguments_.payer.size() && arguments_.legs.size() == arguments_.currency.size(),
               "FmmLsmMultiLegOptionEngine: inconsistent legs / payer / currency vectors");
    for (const auto& c : arguments_.currency)
        QL_REQUIRE(c == p.currency(), "FmmLsmMultiLegOptionEngine: leg currency " << c << " != model currency "
                                                                                  << p.currency());
    QL_REQUIRE(arguments_.settlementType == Settlement::Physical,
               "FmmLsmMultiLegOptionEngine: physical settlement only (cash settlement is not supported)");
    QL_REQUIRE(!arguments_.midCouponExercise, "FmmLsmMultiLegOptionEngine: mid-coupon exercise is not supported");
    QL_REQUIRE(arguments_.noticePeriod == 0 * Days,
               "FmmLsmMultiLegOptionEngine: notice periods are not supported for swaptions in this release");

    // flows from the OPTION HOLDER's view: receive legs with payer = false, pay legs with payer = true
    FmmCallableInstrument inst;
    inst.style = FmmCallableInstrument::Style::Enter;
    FmmLegMapper mapper(*grid_, config_.gridToleranceDays, today);
    for (Size i = 0; i < arguments_.legs.size(); ++i) {
        std::ostringstream what;
        what << "leg " << i;
        mapper.addLeg(arguments_.legs[i], arguments_.payer[i] ? -1.0 : 1.0, inst, what.str());
    }
    QL_REQUIRE(inst.lastFlowIdx >= 1, "FmmLsmMultiLegOptionEngine: no future cash flows");
    const Real underlyingNpv = FmmLegMapper::underlyingCurveValue(inst, p);
    results_.underlyingNpv = underlyingNpv;
    results_.additionalResults["underlyingNpv"] = underlyingNpv;

    if (!arguments_.exercise) {
        // no exercise: the instrument is the underlying itself
        results_.value = underlyingNpv;
        results_.additionalResults["fmmGridMaxDateMismatchDays"] = static_cast<Real>(grid_->maxMismatchDays());
        return;
    }
    QL_REQUIRE(arguments_.exercise->type() != Exercise::American,
               "FmmLsmMultiLegOptionEngine: American exercise is not supported (grid-date Bermudan only)");
    // exercise dates -> rights; every exercise date must be a grid date at which all legs have an
    // accrual boundary, i.e. no coupon straddles it (checked through the flow mapping: a coupon
    // with accrual start before and pay date after an exercise date would be mis-assigned)
    for (const Date& d : arguments_.exercise->dates()) {
        if (d <= today)
            continue;
        const Size idx = grid_->index(d, config_.gridToleranceDays, "exercise date");
        if (idx >= inst.lastFlowIdx)
            continue; // exercise at or after the last flow has no value
        for (const auto& f : inst.compoundedFloats)
            QL_REQUIRE(!(f.startIdx < idx && f.payIdx > idx),
                       "FmmLsmMultiLegOptionEngine: exercise date " << d << " falls inside a compounded coupon period ("
                                                                    << grid_->dates()[f.startIdx] << " - "
                                                                    << grid_->dates()[f.endIdx]
                                                                    << "); whole-period exercise required");
        inst.rights.push_back({idx, idx, 0.0});
    }
    if (inst.rights.empty()) {
        results_.value = 0.0;
        results_.additionalResults["fmmGridMaxDateMismatchDays"] = static_cast<Real>(grid_->maxMismatchDays());
        return;
    }
    // notice indices strictly ascending (duplicates after mapping are collapsed)
    std::vector<FmmCallableInstrument::Right> rights;
    for (const auto& r : inst.rights)
        if (rights.empty() || r.noticeIdx > rights.back().noticeIdx)
            rights.push_back(r);
    inst.rights = rights;

    FmmLsmPricer pricer(model_, inst, config_.lsm);
    const FmmLsmResult res = pricer.calculate();
    results_.value = res.lowerBound;
    std::unique_ptr<FmmDualBoundResult> dual;
    if (config_.dualBound)
        dual = std::make_unique<FmmDualBoundResult>(
            pricer.dualBound(config_.dualOuterPaths, config_.dualInnerPaths, config_.dualSeed));
    fmmWriteLsmResults(results_.additionalResults, res, dual.get(), config_, *grid_);
    results_.additionalResults["fmmExerciseDates"] = arguments_.exercise->dates();
}

} // namespace QuantExt
