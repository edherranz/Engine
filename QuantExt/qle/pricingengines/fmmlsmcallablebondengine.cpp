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

#include <qle/pricingengines/fmmlsmcallablebondengine.hpp>

#include <ql/cashflows/coupon.hpp>

#include <cmath>
#include <cstdlib>
#include <memory>

namespace QuantExt {

FmmLsmCallableBondEngine::FmmLsmCallableBondEngine(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                                                   const QuantLib::ext::shared_ptr<FmmGrid>& grid,
                                                   const FmmLsmEngineConfig& config,
                                                   const Handle<YieldTermStructure>& referenceCurve,
                                                   const Handle<Quote>& discountingSpread,
                                                   const Handle<DefaultProbabilityTermStructure>& creditCurve)
    : model_(model), grid_(grid), config_(config), referenceCurve_(referenceCurve), discountingSpread_(discountingSpread),
      creditCurve_(creditCurve) {
    QL_REQUIRE(model_ && grid_, "FmmLsmCallableBondEngine: model or grid is null");
    QL_REQUIRE(grid_->numberOfRates() == model_->parametrization()->numberOfRates(),
               "FmmLsmCallableBondEngine: grid and model disagree on the number of periods");
    registerWith(model_->parametrization()->termStructure());
    registerWith(referenceCurve_);
    registerWith(discountingSpread_);
    registerWith(creditCurve_);
    for (const auto& o : config_.observables)
        registerWith(o);
}

void FmmLsmCallableBondEngine::calculate() const {
    const auto& p = *model_->parametrization();
    const Date today = grid_->referenceDate();
    QL_REQUIRE(p.termStructure()->referenceDate() == today, "FmmLsmCallableBondEngine: model curve / grid date mismatch");
    QL_REQUIRE(!arguments_.perpetual, "FmmLsmCallableBondEngine: perpetual bonds are not supported");
    QL_REQUIRE(arguments_.putData.empty(), "FmmLsmCallableBondEngine: put rights are not supported in this release");
    QL_REQUIRE(creditCurve_.empty(),
               "FmmLsmCallableBondEngine: credit curves are not supported (deterministic issuer spread only, "
               "FMM_SPEC.md section 8)");
    if (!referenceCurve_.empty())
        QL_REQUIRE(referenceCurve_->referenceDate() == today, "FmmLsmCallableBondEngine: reference curve date mismatch");
    const Real spread = discountingSpread_.empty() ? 0.0 : discountingSpread_->value();

    // issuer's view: every bond cash flow is paid (negative), calls are cancellation rights
    FmmCallableInstrument inst;
    inst.style = FmmCallableInstrument::Style::Cancel;
    inst.issuerSpread = spread;
    FmmLegMapper mapper(*grid_, config_.gridToleranceDays, today);
    mapper.addLeg(arguments_.cashflows, -1.0, inst, "bond cash flows");
    QL_REQUIRE(inst.lastFlowIdx >= 1, "FmmLsmCallableBondEngine: no future cash flows");

    // outstanding notional at a call date: the nominal of the coupon accruing after it (or the
    // last known notional); calls must fall on coupon pay dates, where the accrual is zero
    std::vector<std::pair<Date, Real>> couponDates; // pay date, nominal of that coupon
    for (const auto& cf : arguments_.cashflows)
        if (auto c = QuantLib::ext::dynamic_pointer_cast<Coupon>(cf))
            couponDates.emplace_back(c->date(), c->nominal());
    auto notionalAfter = [&](const Date& d) {
        Real n = arguments_.notionals.empty() ? Null<Real>() : arguments_.notionals.front();
        for (const auto& cf : arguments_.cashflows)
            if (auto c = QuantLib::ext::dynamic_pointer_cast<Coupon>(cf))
                if (c->accrualStartDate() >= d) {
                    n = c->nominal();
                    break;
                }
        QL_REQUIRE(n != Null<Real>(), "FmmLsmCallableBondEngine: cannot determine the notional outstanding at " << d);
        return n;
    };
    for (const auto& cd : arguments_.callData) {
        if (cd.exerciseDate <= today)
            continue;
        QL_REQUIRE(cd.exerciseType == CallableBond::CallabilityData::ExerciseType::OnThisDate,
                   "FmmLsmCallableBondEngine: American (FromThisDateOn) calls are not supported");
        bool onCouponDate = false;
        for (const auto& c : couponDates)
            if (std::abs(static_cast<Integer>(c.first - cd.exerciseDate)) <= static_cast<Integer>(config_.gridToleranceDays))
                onCouponDate = true;
        QL_REQUIRE(onCouponDate, "FmmLsmCallableBondEngine: call date "
                                     << cd.exerciseDate
                                     << " is not a coupon pay date; accrued-interest handling for mid-period "
                                        "calls is not supported in this release");
        const Size idx = grid_->index(cd.exerciseDate, config_.gridToleranceDays, "call date");
        if (idx >= inst.lastFlowIdx)
            continue;
        // clean = dirty on a coupon date (zero accrual); price is per unit notional
        const Real fee = -cd.price * notionalAfter(cd.exerciseDate);
        inst.rights.push_back({idx, idx, fee});
    }
    std::vector<FmmCallableInstrument::Right> rights;
    for (const auto& r : inst.rights)
        if (rights.empty() || r.noticeIdx > rights.back().noticeIdx)
            rights.push_back(r);
    inst.rights = rights;

    const Real issuerStraight = FmmLegMapper::underlyingCurveValue(inst, p); // negative
    Real issuerValue = issuerStraight;
    std::unique_ptr<FmmLsmResult> res;
    std::unique_ptr<FmmDualBoundResult> dual;
    if (!inst.rights.empty()) {
        FmmLsmPricer pricer(model_, inst, config_.lsm);
        res = std::make_unique<FmmLsmResult>(pricer.calculate());
        issuerValue = res->lowerBoundCv; // control variate on the known straight-bond value
        if (config_.dualBound)
            dual = std::make_unique<FmmDualBoundResult>(
                pricer.dualBound(config_.dualOuterPaths, config_.dualInnerPaths, config_.dualSeed));
    }
    // holder's value; the settlement value follows ORE's convention (npv date = today)
    results_.value = -issuerValue;
    const Real settleDf = referenceCurve_.empty() ? 1.0 : referenceCurve_->discount(arguments_.settlementDate);
    results_.settlementValue = results_.value / settleDf;
    results_.additionalResults["strippedBondNpv"] = -issuerStraight;
    results_.additionalResults["callPutValue"] = -issuerStraight - results_.value;
    results_.additionalResults["fmmIssuerSpread"] = spread;
    if (res) {
        fmmWriteLsmResults(results_.additionalResults, *res, dual.get(), config_, *grid_);
        // for the Cancel-style issuer instrument the reported bounds are the control-variate ones
        results_.additionalResults["fmmLsmLowerBound"] = res->lowerBoundCv;
        results_.additionalResults["fmmLsmLowerBoundStdError"] = res->lowerBoundCvSe;
        if (dual) {
            results_.additionalResults["fmmLsmUpperBound"] = res->lowerBoundCv + dual->gap;
            results_.additionalResults["fmmLsmUpperBoundStdError"] =
                std::sqrt(res->lowerBoundCvSe * res->lowerBoundCvSe + dual->gapSe * dual->gapSe);
        }
        results_.additionalResults["fmmValuationPerspective"] =
            std::string("issuer (Cancel style); NPV reported for the holder as -(issuer value)");
    }
    if (config_.calibrationResults)
        for (const auto& kv : config_.calibrationResults())
            results_.additionalResults[kv.first] = kv.second;
}

} // namespace QuantExt
