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

/*! \file fmmgrid.hpp
    \brief FMM tenor grid on calendar dates: construction from contractual dates, date-to-index
           mapping with an explicit tolerance (no silent snapping), swap specs from QuantLib
           swaption arguments (ore-fmm docs/A5_PLAN.md section 3)
    \ingroup models
*/

#ifndef quantext_fmm_grid_hpp
#define quantext_fmm_grid_hpp

#include <qle/models/fmmanalytics.hpp>

#include <ql/instruments/swaption.hpp>
#include <ql/time/businessdayconvention.hpp>
#include <ql/time/calendar.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/date.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/period.hpp>

#include <set>
#include <string>
#include <vector>

namespace QuantExt {
using namespace QuantLib;

//! tenor grid T_0 = reference date < T_1 < ... < T_M on dates, with model times in the day counter
class FmmGrid {
public:
    FmmGrid() = default;
    //! grid = reference date, the tenor lattice anchored on the last contractual date (walking
    //! backwards to the reference date), the contractual dates (sorted, unique, strictly after the
    //! reference date) and, with a notice period, the notice date of every grid date. Lattice and
    //! notice dates within mergeToleranceDays of a contractual date are dropped in its favour, so
    //! that no two generated grid dates lie within the engines' date-mapping tolerance.
    FmmGrid(const Date& referenceDate, const std::set<Date>& contractualDates, const Period& fillerTenor,
            const DayCounter& dc = Actual365Fixed(), const Natural minimumStubDays = 2,
            const Period& noticePeriod = Period(0, Days), const Calendar& noticeCalendar = NullCalendar(),
            const BusinessDayConvention noticeConvention = Preceding, const Natural mergeToleranceDays = 3);

    const Date& referenceDate() const { return referenceDate_; }
    const std::vector<Date>& dates() const { return dates_; }
    const Array& times() const { return times_; }
    Size numberOfRates() const { return dates_.size() - 1; }
    const DayCounter& dayCounter() const { return dc_; }

    //! index of the grid date nearest to d, provided it lies within toleranceDays; otherwise a
    //! descriptive exception naming `what`. The largest accepted mismatch is tracked for reporting.
    Size index(const Date& d, const Natural toleranceDays, const std::string& what) const;
    //! largest date mismatch (in days) accepted by index() so far
    Natural maxMismatchDays() const { return maxMismatchDays_; }

private:
    Date referenceDate_;
    std::vector<Date> dates_;
    Array times_;
    DayCounter dc_;
    mutable Natural maxMismatchDays_ = 0;
};

//! swap spec of the underlying of a physically settled swaption on the grid: expiry = exercise
//! date, float leg = grid periods (a, b], fixed leg from the swap's fixed pay dates and coupon
//! accruals (recovered from the coupon amounts, nominal and rate). Rejects mid-period exercise
//! and fixed/float dates off the grid beyond the tolerance. With a model curve, the floating
//! leg's time-0 forward basis over that curve (a term-rate index such as EURIBOR, or a payment
//! convention, against the grid's par RFR leg) is recorded as deterministic basis amounts, so that
//! the model's forward swap rate equals the helper's market forward discounted on the model curve.
FmmSwapSpec fmmSwapSpecFromSwaption(const FmmGrid& grid, const Swaption::arguments& args, const Natural toleranceDays,
                                    Real& strike, Option::Type& type, Real& nominal,
                                    const Handle<YieldTermStructure>& modelCurve = Handle<YieldTermStructure>());

} // namespace QuantExt

#endif
