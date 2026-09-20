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

/*! \file portfolio/builders/fmm.hpp
    \brief FMM engine builders (model "FMM", engine "LSM") for swaptions and callable bonds, and the
           per-trade model construction from pricingengine.xml ModelParameters (ore-fmm
           docs/A5_PLAN.md section 6)
    \ingroup builders
*/

#pragma once

#include <ored/model/fmmbuilder.hpp>
#include <ored/portfolio/builders/callablebond.hpp>
#include <ored/portfolio/builders/swaption.hpp>

#include <qle/pricingengines/fmmlsmmultilegoptionengine.hpp>

namespace ore {
namespace data {

//! what a trade needs from the model: calibration basket definition and the dates the grid must hold
struct FmmModelRequest {
    std::string id;
    std::string key; //!< index name or currency (the IrModelBuilder qualifier)
    std::vector<Date> expiries, maturities;
    std::vector<Real> strikes; //!< Null<Real>() = ATM
    std::vector<Date> gridDates;
};

/*! ModelParameters (all under <ModelParameters>): Calibration, CalibrationStrategy, Tolerance,
    VolatilityType, Shift, RhoInf, Beta, Factors, Volatility, VolatilityTimes, ApproximationMethod,
    McCorrected, McPathsPerRep, McReps, McSeed, McMaxIterations, McToleranceBp, Grid,
    GridToleranceDays, SubSteps. The builder is registered with the engine factory's model builders
    so that scenario runs recalibrate it through the same lazy path as the LGM builders. */
QuantLib::ext::shared_ptr<FmmBuilder> fmmModelBuilder(const EngineBuilder* builder, const FmmModelRequest& request);

/*! EngineParameters: TrainingPaths, ValuationPaths, TrainingSeed, ValuationSeed, TrainingSequence,
    ValuationSequence, BasisOrder, DualBound, DualOuterPaths, DualInnerPaths, DualSeed */
QuantExt::FmmLsmEngineConfig fmmLsmEngineConfig(const EngineBuilder* builder,
                                                const QuantLib::ext::shared_ptr<FmmBuilder>& fmmBuilder);

//! Bermudan / European swaptions on the FMM LSM engine
class FmmLsmSwaptionEngineBuilder final : public SwaptionEngineBuilder {
public:
    FmmLsmSwaptionEngineBuilder()
        : SwaptionEngineBuilder("FMM", "LSM",
                                {"EuropeanSwaption", "EuropeanSwaption_NonStandard", "BermudanSwaption",
                                 "BermudanSwaption_NonStandard"}) {}

private:
    QuantLib::ext::shared_ptr<PricingEngine>
    engineImpl(const string& id, const std::vector<string>& keys, const std::vector<Date>& dates,
               const std::vector<Date>& maturities, const std::vector<std::vector<Real>>& strikes,
               const std::vector<std::vector<Real>>& fxStrikes, const bool isAmerican, const std::string& discountCurve,
               const std::string& securitySpread, const CamOrLgmModel&) override;
};

//! callable bonds on the FMM LSM engine (deterministic issuer spread from the security spread)
class FmmLsmCallableBondEngineBuilder final : public CallableBondEngineBuilder {
public:
    FmmLsmCallableBondEngineBuilder() : CallableBondEngineBuilder("FMM", "LSM") {}

protected:
    QuantLib::ext::shared_ptr<QuantExt::PricingEngine>
    engineImpl(const std::string& id, const std::string& ccy, const std::string& creditCurveId,
               const std::string& securityId, const std::string& referenceCurveId, const std::string& incomeCurveId,
               const QuantLib::Date& maturityDate) override;
};

} // namespace data
} // namespace ore
