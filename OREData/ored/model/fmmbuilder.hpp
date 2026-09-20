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

/*! \file ored/model/fmmbuilder.hpp
    \brief Build and calibrate a generalized Forward Market Model (ore-fmm docs/A5_PLAN.md
           section 6): inherits the swaption basket, market observation and lazy recalibration
           of IrModelBuilder; the model itself is standalone (no CrossAssetModel integration)
    \ingroup models
*/

#pragma once

#include <ored/marketdata/market.hpp>
#include <ored/model/fmmdata.hpp>
#include <ored/model/irmodelbuilder.hpp>

#include <ql/quotes/simplequote.hpp>

#include <qle/models/fmmgrid.hpp>
#include <qle/models/fmmirmodel.hpp>
#include <qle/models/forwardmarketmodel.hpp>

#include <ql/any.hpp>

#include <map>

namespace ore {
namespace data {

//! calibration record (model-neutral fields + FMM specifics), exported as additional results
struct FmmCalibrationInfo {
    bool valid = false;
    Real rmse = Null<Real>();                //!< IrModelBuilder's calibration error over the basket
    std::string approximationMethod;
    std::vector<QuantExt::SwaptionData> swaptionData; //!< basket details (inherited format)
    std::vector<Real> expiryTimes, strikes, marketVols, modelVols, effectiveShifts, mcVols, mcVolStdErrors;
    //! cap/floor basket (joint bootstrap), empty when no basket was configured
    std::vector<Real> capletTimes, capletStrikes, capletMarketVols, capletModelVols;
    Size jointIterations = 0;
    bool jointConverged = false;  //!< both baskets within the joint tolerance
    bool jointStationary = false; //!< the alternation reached its fixed point
    Real jointLastParameterChange = Null<Real>();
    Real capletMaxResidualBp = Null<Real>();
    Size gridPeriods = 0;
    Natural gridMaxDateMismatchDays = 0;
    Size mcIterations = 0;
    bool mcConverged = false;
    std::map<std::string, QuantLib::ext::any> additionalResults() const;
};

class FmmBuilder : public IrModelBuilder {
public:
    FmmBuilder(const QuantLib::ext::shared_ptr<ore::data::Market>& market, const QuantLib::ext::shared_ptr<FmmData>& data,
               const std::string& configuration = Market::defaultConfiguration, Real bootstrapTolerance = 0.001,
               const bool continueOnError = false, const std::string& referenceCalibrationGrid = "",
               const bool setCalibrationInfo = false, const std::string& id = "unknown",
               BlackCalibrationHelper::CalibrationErrorType calibrationErrorType =
                   BlackCalibrationHelper::RelativePriceError,
               const bool allowChangingFallbacksUnderScenarios = false, const bool allowModelFallbacks = false,
               const bool dontCalibrate = false);

    //! the calibrated model (triggers the lazy calculation)
    QuantLib::ext::shared_ptr<QuantExt::ForwardMarketModel> modelAsFmm() const;
    //! the model's date grid (triggers the lazy calculation)
    QuantLib::ext::shared_ptr<QuantExt::FmmGrid> grid() const;
    //! the IrModel adapter (observable: notifies after each recalibration)
    QuantLib::ext::shared_ptr<QuantExt::FmmIrModel> irModel() const;
    const FmmCalibrationInfo& calibrationInfo() const { return info_; }

private:
    void initParametrization() const override;
    void calibrate() const override;
    QuantLib::ext::shared_ptr<PricingEngine> getPricingEngine() const override;

    bool setCalibrationInfo_ = false;
    mutable Date parametrizationInitializedOnAnchorDate_;
    mutable QuantLib::ext::shared_ptr<QuantExt::FmmGrid> grid_;
    mutable std::vector<Real> levelScale_; //!< tau_j / tau_grid per period (level volatility scaling)
    mutable QuantLib::ext::shared_ptr<QuantExt::FmmParametrization> fmmParametrization_;
    mutable QuantLib::ext::shared_ptr<QuantExt::ForwardMarketModel> fmm_;
    mutable QuantLib::ext::shared_ptr<QuantExt::FmmIrModel> irModel_;
    mutable FmmCalibrationInfo info_;
};

} // namespace data
} // namespace ore
