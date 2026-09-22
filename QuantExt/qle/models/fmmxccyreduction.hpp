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

/*! \file qle/models/fmmxccyreduction.hpp
    \brief Hagan's reduction of a cancellable cross-currency swap to a single-currency callable
           structure on the FMM (Hagan, Callable X-CCY Swap, sections 1.1, 2, 3.1-3.2; ore-fmm
           docs/A9_1_HAGAN_REDUCTION_PLAN.md, A9_XCCY_TEMPLATE_PLAN.md)

    The deal is written as an ORE MultiLegOption (a Bermudan option into a two-currency swap with
    exercise rebates). Currency A is the floating currency, telescoped away by the bank-account
    identity [1 + tau r] Z(t; T_{i-1}) = Z(t; T_i) (eq. 2.2): the rate part of every A floating
    coupon is worth +N at its accrual start and -N at its accrual end, so the A side collapses to
    deterministic amounts (spread stream, notional flows, rebates, coupon starts) that today's
    forward FX F^{A/B}(0; T) = F^{A/B}(0) Z^A(0; T) / Z^B(0; T) (eq. 1.7b) converts into currency
    B. The B side is kept as it is. The result is an Enter-style FmmCallableInstrument in B whose
    exercise payoff is N_B L_k^B [R_fix - s_k^eq - R_k^B] (eq. 1.8a) with the equivalent spread
    s_k^eq (eq. 1.8d / 3.3b); the per-right diagnostics carry the calibration strikes
    R_fix - s_k^eq (eq. 3.3a), the annuity / convexity ratios gamma (eqs. 3.12, 3.15) and the
    quantities the improved method's effective volatility (eq. 3.41) needs.
*/

#pragma once

#include <qle/instruments/multilegoption.hpp>
#include <qle/models/fmmgrid.hpp>
#include <qle/models/fmmlsmpricer.hpp>

#include <ql/currency.hpp>
#include <ql/handle.hpp>
#include <ql/quote.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>

#include <set>
#include <vector>

namespace QuantExt {
using namespace QuantLib;

//! curves and spot of the reduction (collateral-consistent pricing curves, Hagan (f))
struct FmmXccyReductionInputs {
    Currency ccyA; //!< the telescoped currency (Hagan's ccy A: the floating leg with the spread)
    Currency ccyB; //!< the modelled currency (Hagan's ccy B: the fixed leg)
    Handle<YieldTermStructure> discountA, discountB; //!< pricing discount curves of A and B
    Handle<Quote> fxSpotBperA;                       //!< units of B per unit of A (Hagan's F^{A/B}(t))
};

//! per exercise right: the reduced B-currency swaption geometry (Hagan section 3)
struct FmmXccyReductionRight {
    Date noticeDate, settleDate;
    Size noticeIdx = 0, settleIdx = 0;
    Time expiry = 0.0;                         //!< notice time t_k^ex (years)
    Time swapLengthA = 0.0, swapLengthB = 0.0; //!< settle -> last A / B payment (years)
    Real notionalA = 0.0, notionalB = 0.0;
    Real annuityA = 0.0, annuityB = 0.0; //!< L_k^A(0), L_k^B(0) per unit notional, own accruals (eqs. 1.8c, 1.8e)
    Real gammaA = 0.0, gammaB = 0.0;     //!< C^{fwd} / L^{fwd} from the notice date (eqs. 3.12b, 3.15b)
    Real forwardSwapRateA = 0.0;         //!< ATM rate of the A coterminal swap on the A discount curve
    Real fixedRate = 0.0;                //!< annuity-weighted B fixed rate after settlement (R_fix^B)
    Real equivalentSpread = 0.0;         //!< s_k^eq (eq. 3.3b), positive when the spread offsets the fixed rate
    Real strike = 0.0;                   //!< R_fix - s_k^eq (eq. 3.3a)
    Real forwardFxAtSettle = 0.0;        //!< F^{A/B}(0; T_settle)
    Real feeFlowB = 0.0;                 //!< the right's fee in B (rebates and the converted A coupon start)
};

struct FmmXccyReductionDiagnostics {
    std::vector<FmmXccyReductionRight> rights;
    Real fxSpotBperA = 0.0;
    Real spreadStreamValueB = 0.0; //!< t = 0 B value of the converted A spread stream (all future coupons)
    Size telescopedCoupons = 0;    //!< A floating coupons telescoped
    Size convertedFlows = 0;       //!< deterministic A flows converted with the forward FX
};

//! contractual dates of the legs and the exercise (accrual boundaries, payments, notice and
//! settlement dates after today): the model grid must carry them
std::set<Date> fmmXccyContractualDates(const MultiLegOption::arguments& args, const Date& today);

/*! reduced instrument in currency B from the option holder's view. Eligibility (hard rejections):
    every A floating coupon is a compounded RFR coupon whose index projects on the A discount
    curve (the telescoping identity), with unit gearing, no compounded spread, no lookback /
    cut-off / lockout, payment at the accrual end; whole periods at every settlement date (no
    coupon of either leg straddles it); rebates paid on the exercise date, in A or B. */
FmmCallableInstrument fmmMapXccyReduced(const MultiLegOption::arguments& args, const FmmGrid& grid,
                                        const Natural toleranceDays, const Date& today,
                                        const FmmXccyReductionInputs& inputs, FmmXccyReductionDiagnostics& diag,
                                        std::vector<Date>& usedSettle);

} // namespace QuantExt
