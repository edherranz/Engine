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

#include <qle/models/fmmgrid.hpp>

#include <ql/errors.hpp>
#include <ql/exercise.hpp>
#include <ql/instruments/fixedvsfloatingswap.hpp>

#include <algorithm>
#include <cstdlib>
#include <iterator>

namespace QuantExt {

FmmGrid::FmmGrid(const Date& referenceDate, const std::set<Date>& contractualDates, const Period& fillerTenor,
                 const DayCounter& dc, const Natural minimumStubDays)
    : referenceDate_(referenceDate), dc_(dc) {
    QL_REQUIRE(fillerTenor.length() > 0, "FmmGrid: filler tenor must be positive");
    std::set<Date> all;
    for (const Date& d : contractualDates)
        if (d > referenceDate)
            all.insert(d);
    QL_REQUIRE(!all.empty(), "FmmGrid: no contractual date after the reference date " << referenceDate);
    // fill backwards from the first contractual date with the filler tenor (the front stub at
    // the reference date is at most one tenor, and no grid date is closer than minimumStubDays
    // to it), then fill every gap between consecutive contractual dates with tenor steps from
    // the earlier date, so that the model grid keeps the tenor's resolution (e.g. quarterly rates
    // under annual coupons) while every contractual date stays an exact grid date
    const Date first = *all.begin();
    for (Size k = 1; k < 10000; ++k) {
        const Date f = first - static_cast<Integer>(k) * fillerTenor;
        if (f <= referenceDate + static_cast<Integer>(minimumStubDays))
            break;
        all.insert(f);
    }
    const std::vector<Date> contractual(all.begin(), all.end());
    for (Size i = 0; i + 1 < contractual.size(); ++i) {
        for (Size k = 1; k < 10000; ++k) {
            const Date f = contractual[i] + static_cast<Integer>(k) * fillerTenor;
            if (f >= contractual[i + 1] - static_cast<Integer>(minimumStubDays))
                break;
            all.insert(f);
        }
    }
    dates_.push_back(referenceDate);
    for (const Date& x : all)
        if (x - referenceDate >= static_cast<Integer>(minimumStubDays))
            dates_.push_back(x);
    QL_REQUIRE(dates_.size() >= 2, "FmmGrid: grid needs at least one period");
    times_ = Array(dates_.size());
    for (Size k = 0; k < dates_.size(); ++k)
        times_[k] = dc_.yearFraction(referenceDate_, dates_[k]);
    for (Size k = 1; k < dates_.size(); ++k)
        QL_REQUIRE(times_[k] > times_[k - 1], "FmmGrid: non-increasing model times at " << dates_[k]);
}

Size FmmGrid::index(const Date& d, const Natural toleranceDays, const std::string& what) const {
    auto it = std::lower_bound(dates_.begin(), dates_.end(), d);
    Size best = Null<Size>();
    Integer bestDiff = 0;
    auto consider = [&](std::vector<Date>::const_iterator c) {
        if (c == dates_.end())
            return;
        const Integer diff = std::abs(static_cast<Integer>(*c - d));
        if (best == Null<Size>() || diff < bestDiff) {
            best = static_cast<Size>(c - dates_.begin());
            bestDiff = diff;
        }
    };
    consider(it);
    if (it != dates_.begin())
        consider(std::prev(it));
    QL_REQUIRE(best != Null<Size>(), "FmmGrid: empty grid");
    QL_REQUIRE(bestDiff <= static_cast<Integer>(toleranceDays),
               "FmmGrid: " << what << " " << d << " is " << bestDiff << " days from the nearest grid date "
                           << dates_[best] << " (tolerance " << toleranceDays
                           << " days) - the structure is not representable on the model grid");
    maxMismatchDays_ = std::max(maxMismatchDays_, static_cast<Natural>(bestDiff));
    return best;
}

FmmSwapSpec fmmSwapSpecFromSwaption(const FmmGrid& grid, const Swaption::arguments& args, const Natural toleranceDays,
                                    Real& strike, Option::Type& type, Real& nominal) {
    QL_REQUIRE(args.swap, "fmmSwapSpecFromSwaption: swaption without underlying swap");
    QL_REQUIRE(args.exercise && args.exercise->dates().size() == 1,
               "fmmSwapSpecFromSwaption: a single (European) exercise date is required");
    const auto& swap = *args.swap;
    nominal = swap.nominal();
    strike = swap.fixedRate();
    type = swap.type() == Swap::Payer ? Option::Call : Option::Put;
    const std::vector<Date>& floatDates = swap.floatingSchedule().dates();
    QL_REQUIRE(floatDates.size() >= 2, "fmmSwapSpecFromSwaption: degenerate floating schedule");
    FmmSwapSpec spec;
    // the swap starts at the first floating accrual date; the exercise date must coincide with it
    // (grid-aligned expiry) within the tolerance
    spec.a = grid.index(floatDates.front(), toleranceDays, "swap start date");
    const Size exIdx = grid.index(args.exercise->dates().front(), toleranceDays, "exercise date");
    QL_REQUIRE(exIdx == spec.a, "fmmSwapSpecFromSwaption: exercise date " << args.exercise->dates().front()
                                                                           << " does not coincide with the swap start "
                                                                           << floatDates.front()
                                                                           << " (mid-period exercise is unsupported)");
    spec.b = grid.index(floatDates.back(), toleranceDays, "swap end date");
    QL_REQUIRE(spec.b > spec.a, "fmmSwapSpecFromSwaption: swap end not after start on the grid");
    // fixed leg: pay dates and accruals recovered from coupon amounts (nominal * rate * accrual)
    const std::vector<Date>& fixedDates = swap.fixedSchedule().dates();
    const DayCounter& fdc = swap.fixedDayCount();
    for (Size c = 1; c < fixedDates.size(); ++c) {
        spec.fixedPayIndices.push_back(grid.index(fixedDates[c], toleranceDays, "fixed pay date"));
        spec.fixedAccruals.push_back(fdc.yearFraction(fixedDates[c - 1], fixedDates[c]));
    }
    spec.validate(grid.numberOfRates());
    return spec;
}

} // namespace QuantExt
