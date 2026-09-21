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

/*! \file ored/model/fmmdata.hpp
    \brief Generalized Forward Market Model data (ore-fmm docs/A5_PLAN.md section 6): model
           parameters, calibration instructions and the approximation / MC-correction switches
    \ingroup models
*/

#pragma once

#include <ored/marketdata/market.hpp>
#include <ored/model/irmodeldata.hpp>

#include <qle/models/fmmanalytics.hpp>
#include <qle/models/fmmparametrization.hpp>

#include <ql/time/businessdayconvention.hpp>
#include <ql/time/calendar.hpp>
#include <ql/time/date.hpp>
#include <ql/time/period.hpp>

#include <vector>

namespace ore {
namespace data {
using namespace QuantLib;

class FmmData : public IrModelData {
public:
    struct McCorrection {
        bool enabled = false;
        Size pathsPerRep = 16384;
        Size reps = 8;
        BigNatural seed = 9100;
        Size maxIterations = 3;
        Real toleranceBp = 0.2;
    };

    FmmData();

    void clear() override;
    void reset() override;

    //! \name Serialisation
    //@{
    void fromXML(XMLNode* node) override;
    XMLNode* toXML(XMLDocument& doc) const override;
    //@}

    //! consistency checks (throws)
    void validate() const;

    //! \name Setters/Getters
    //@{
    Period& grid() { return grid_; }
    //! displaced-diffusion shift; Null<Real>() means 1/tau_j (the HJM replication shift)
    Real& shift() { return shift_; }
    QuantExt::FmmParametrization::LocalVolType& volatilityType() { return volType_; }
    Real& rhoInf() { return rhoInf_; }
    Real& beta() { return beta_; }
    Size& factors() { return factors_; }
    bool& calibrateVolatility() { return calibrateVol_; }
    ParamType& volParamType() { return volParamType_; }
    std::vector<Time>& volTimes() { return volTimes_; }
    std::vector<Real>& volValues() { return volValues_; }
    QuantExt::FmmSwaptionApproxMethod& approximationMethod() { return approx_; }
    McCorrection& mcCorrection() { return mc_; }
    Size& subSteps() { return subSteps_; }
    Natural& gridToleranceDays() { return gridToleranceDays_; }
    //! notice period added to the model grid (the notice date of every grid date), for products whose
    //! notice dates the engine builder cannot see (callable bonds); 0D = none
    Period& noticePeriod() { return noticePeriod_; }
    Calendar& noticeCalendar() { return noticeCalendar_; }
    BusinessDayConvention& noticeConvention() { return noticeConvention_; }
    //! cap/floor (optionlet) calibration basket: "None" or "ATM" (one backward-looking caplet per
    //! grid period at the period's forward, target vol from the market's stripped optionlet surface
    //! keyed by the qualifier), fitted jointly with the swaption basket (A3 joint bootstrap)
    std::string& capFloorBasket() { return capFloorBasket_; }
    //! last period end included in the cap/floor basket (0D = the whole grid)
    Period& capFloorHorizon() { return capFloorHorizon_; }
    Size& jointMaxIterations() { return jointMaxIterations_; }
    Real& jointToleranceBp() { return jointToleranceBp_; }
    //! CalibrationType BestFit only: fit the a(t) segments (one per distinct swaption expiry) by
    //! least squares to the whole swaption basket, e.g. an expiry x term grid, jointly with the
    //! cap/floor basket when one is configured (A7 scenario model). Default false = the common
    //! level fit (mean residual zero across the basket, a(t) as given).
    bool& bestFitTimeDependence() { return bestFitTimeDependence_; }
    //! contractual dates the model grid must contain (set by engine builders from the trade;
    //! not part of the XML). Empty: the grid is built from the calibration basket dates.
    std::vector<Date>& gridDates() { return gridDates_; }
    std::vector<std::string>& optionExpiries() const { return optionExpiries_; }
    std::vector<std::string>& optionTerms() const { return optionTerms_; }
    std::vector<std::string>& optionStrikes() const { return optionStrikes_; }
    //@}

private:
    Period grid_;
    Real shift_;
    QuantExt::FmmParametrization::LocalVolType volType_;
    Real rhoInf_, beta_;
    Size factors_;
    bool calibrateVol_;
    ParamType volParamType_;
    std::vector<Time> volTimes_;
    std::vector<Real> volValues_;
    QuantExt::FmmSwaptionApproxMethod approx_;
    McCorrection mc_;
    Size subSteps_;
    Natural gridToleranceDays_;
    Period noticePeriod_;
    Calendar noticeCalendar_;
    BusinessDayConvention noticeConvention_;
    std::string capFloorBasket_;
    Period capFloorHorizon_;
    Size jointMaxIterations_;
    Real jointToleranceBp_;
    bool bestFitTimeDependence_ = false;
    std::vector<Date> gridDates_;
    mutable std::vector<std::string> optionExpiries_, optionTerms_, optionStrikes_;
};

QuantExt::FmmParametrization::LocalVolType parseFmmVolatilityType(const std::string& s);
QuantExt::FmmSwaptionApproxMethod parseFmmApproximationMethod(const std::string& s);
//! "1/tau" or a number
Real parseFmmShift(const std::string& s);
std::string to_string(const QuantExt::FmmParametrization::LocalVolType& t);
std::string to_string(const QuantExt::FmmSwaptionApproxMethod& m);

} // namespace data
} // namespace ore
