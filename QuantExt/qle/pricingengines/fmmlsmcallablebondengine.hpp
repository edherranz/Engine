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

/*! \file fmmlsmcallablebondengine.hpp
    \brief CallableBond engine on the FMM Longstaff-Schwartz pricer: issuer calls on coupon dates
           (clean = dirty there), deterministic issuer spread, holder NPV = -(issuer value), with
           the stripped bond and embedded option as additional results (ore-fmm docs/A5_PLAN.md
           section 5). Puts, American calls, credit curves and off-grid call dates are rejected.
    \ingroup engines
*/

#ifndef quantext_fmm_lsm_callable_bond_engine_hpp
#define quantext_fmm_lsm_callable_bond_engine_hpp

#include <qle/instruments/callablebond.hpp>
#include <qle/pricingengines/fmmlsmmultilegoptionengine.hpp>

#include <ql/optional.hpp>

#include <ql/quote.hpp>
#include <ql/termstructures/defaulttermstructure.hpp>

namespace QuantExt {

//! CallableBond arguments onto the ISSUER's Cancel-style instrument on the grid: every bond cash
//! flow paid (negative), each future call a cancellation right (decision at the notice date, fee =
//! -price x outstanding notional at the call date, calls on coupon pay dates only), the deterministic
//! issuer spread as the discounting spread. noticeDates receives the decision dates of the rights.
FmmCallableInstrument fmmMapCallableBond(const CallableBond::arguments& args, const FmmGrid& grid,
                                         const Natural toleranceDays, const Date& today, const Real issuerSpread,
                                         std::vector<Date>& noticeDates);

class FmmLsmCallableBondEngine : public CallableBond::engine {
public:
    FmmLsmCallableBondEngine(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                             const QuantLib::ext::shared_ptr<FmmGrid>& grid, const FmmLsmEngineConfig& config,
                             const Handle<YieldTermStructure>& referenceCurve,
                             const Handle<Quote>& discountingSpread = Handle<Quote>(),
                             const Handle<DefaultProbabilityTermStructure>& creditCurve =
                                 Handle<DefaultProbabilityTermStructure>());
    void calculate() const override;

private:
    QuantLib::ext::shared_ptr<ForwardMarketModel> model_;
    QuantLib::ext::shared_ptr<FmmGrid> grid_;
    FmmLsmEngineConfig config_;
    mutable QuantLib::ext::optional<FmmLsmPolicy> frozenPolicy_;
    Handle<YieldTermStructure> referenceCurve_;
    Handle<Quote> discountingSpread_;
    Handle<DefaultProbabilityTermStructure> creditCurve_;
};

} // namespace QuantExt

#endif
