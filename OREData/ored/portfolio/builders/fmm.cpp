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

#include <ored/portfolio/builders/fmm.hpp>

#include <ql/time/calendars/nullcalendar.hpp>
#include <ored/utilities/dategrid.hpp>
#include <ored/utilities/log.hpp>
#include <ored/utilities/parsers.hpp>
#include <ored/utilities/to_string.hpp>

#include <qle/pricingengines/fmmlsmcallablebondengine.hpp>

#include <algorithm>
#include <set>

using namespace QuantLib;
using namespace QuantExt;

namespace ore {
namespace data {

QuantLib::ext::shared_ptr<FmmBuilder> fmmModelBuilder(const EngineBuilder* builder, const FmmModelRequest& req) {
    QL_REQUIRE(req.expiries.size() == req.maturities.size() && req.strikes.size() == req.expiries.size(),
               "fmmModelBuilder(): expiries (" << req.expiries.size() << "), maturities (" << req.maturities.size()
                                              << "), strikes (" << req.strikes.size() << ") mismatch");

    const auto calibration = parseCalibrationType(builder->modelParameter("Calibration"));
    const auto calibrationStrategy = parseCalibrationStrategy(builder->modelParameter("CalibrationStrategy"));
    const std::vector<std::pair<CalibrationType, CalibrationStrategy>> validCalPairs = {
        {CalibrationType::None, CalibrationStrategy::None},
        {CalibrationType::Bootstrap, CalibrationStrategy::CoterminalATM},
        {CalibrationType::Bootstrap, CalibrationStrategy::CoterminalDealStrike},
        {CalibrationType::BestFit, CalibrationStrategy::CoterminalATM},
        {CalibrationType::BestFit, CalibrationStrategy::CoterminalDealStrike}};
    QL_REQUIRE(std::find(validCalPairs.begin(), validCalPairs.end(), std::make_pair(calibration, calibrationStrategy)) !=
                   validCalPairs.end(),
               "FMM: Calibration (" << calibration << ") and CalibrationStrategy (" << calibrationStrategy
                                    << ") are not allowed in this combination");

    auto data = QuantLib::ext::make_shared<FmmData>();
    data->reset();
    data->qualifier() = req.key;
    data->calibrationType() = calibration;
    data->calibrateVolatility() = calibration != CalibrationType::None;
    data->grid() = parsePeriod(builder->modelParameter("Grid", {}, false, "3M"));
    data->shift() = parseFmmShift(builder->modelParameter("Shift", {}, false, "1/tau"));
    data->volatilityType() =
        parseFmmVolatilityType(builder->modelParameter("VolatilityType", {}, false, "DisplacedDiffusion"));
    data->rhoInf() = parseReal(builder->modelParameter("RhoInf", {}, false, "0.6"));
    data->beta() = parseReal(builder->modelParameter("Beta", {}, false, "0.08"));
    data->factors() = static_cast<Size>(parseInteger(builder->modelParameter("Factors", {}, false, "3")));
    data->volParamType() = ParamType::Piecewise;
    data->volValues() = parseListOfValues<Real>(builder->modelParameter("Volatility", {req.key}), &parseReal);
    data->volTimes() =
        parseListOfValues<Real>(builder->modelParameter("VolatilityTimes", {req.key}, false), &parseReal);
    data->approximationMethod() =
        parseFmmApproximationMethod(builder->modelParameter("ApproximationMethod", {}, false, "EffectiveShift"));
    data->mcCorrection().enabled = parseBool(builder->modelParameter("McCorrected", {}, false, "false"));
    data->mcCorrection().pathsPerRep =
        static_cast<Size>(parseInteger(builder->modelParameter("McPathsPerRep", {}, false, "16384")));
    data->mcCorrection().reps = static_cast<Size>(parseInteger(builder->modelParameter("McReps", {}, false, "8")));
    data->mcCorrection().seed = static_cast<BigNatural>(parseInteger(builder->modelParameter("McSeed", {}, false, "9100")));
    data->mcCorrection().maxIterations =
        static_cast<Size>(parseInteger(builder->modelParameter("McMaxIterations", {}, false, "3")));
    data->mcCorrection().toleranceBp = parseReal(builder->modelParameter("McToleranceBp", {}, false, "0.2"));
    data->subSteps() = static_cast<Size>(parseInteger(builder->modelParameter("SubSteps", {}, false, "1")));
    data->gridToleranceDays() =
        static_cast<Natural>(parseInteger(builder->modelParameter("GridToleranceDays", {}, false, "3")));
    // notice period of the trade's calls (the builder cannot see the trade's call data): the grid
    // then carries the notice date of every grid date, onto which the engine maps the trade's
    data->noticePeriod() = parsePeriod(builder->modelParameter("NoticePeriod", {}, false, "0D"));
    data->noticeCalendar() = parseCalendar(builder->modelParameter("NoticeCalendar", {}, false, "NullCalendar"));
    data->noticeConvention() =
        parseBusinessDayConvention(builder->modelParameter("NoticeConvention", {}, false, "Preceding"));
    // optional cap/floor (optionlet) basket fitted jointly with the swaptions
    data->capFloorBasket() = builder->modelParameter("CapFloorBasket", {}, false, "None");
    data->capFloorHorizon() = parsePeriod(builder->modelParameter("CapFloorHorizon", {}, false, "0D"));
    data->jointMaxIterations() =
        static_cast<Size>(parseInteger(builder->modelParameter("JointMaxIterations", {}, false, "50")));
    data->jointToleranceBp() = parseReal(builder->modelParameter("JointToleranceBp", {}, false, "0.1"));

    if (calibrationStrategy != CalibrationStrategy::None) {
        for (Size i = 0; i < req.expiries.size(); ++i) {
            data->optionExpiries().push_back(ore::data::to_string(req.expiries[i]));
            data->optionTerms().push_back(ore::data::to_string(req.maturities[i]));
            const bool dealStrike = calibrationStrategy == CalibrationStrategy::CoterminalDealStrike &&
                                    req.strikes[i] != Null<Real>();
            data->optionStrikes().push_back(dealStrike ? std::to_string(req.strikes[i]) : std::string("ATM"));
        }
    }
    data->gridDates() = req.gridDates;
    data->validate();

    const Real tolerance = parseReal(builder->modelParameter("Tolerance", {}, false, "0.001"));
    const bool continueOnCalibrationError = builder->globalParameters().count("ContinueOnCalibrationError") > 0 &&
                                            parseBool(builder->globalParameters().at("ContinueOnCalibrationError"));
    const bool allowModelFallbacks = builder->globalParameters().count("AllowModelFallbacks") > 0 &&
                                     parseBool(builder->globalParameters().at("AllowModelFallbacks"));
    auto rt = builder->globalParameters().find("RunType");
    const bool allowChangingFallbacks = rt != builder->globalParameters().end() && rt->second != "SensitivityDelta" &&
                                        rt->second != "SensitivityDeltaGamma";
    auto dc = builder->globalParameters().find("Calibrate");
    const bool dontCalibrate = dc != builder->globalParameters().end() && !parseBool(dc->second);

    auto fmmBuilder = QuantLib::ext::make_shared<FmmBuilder>(
        builder->market(), data, builder->configuration(MarketContext::irCalibration), tolerance,
        continueOnCalibrationError, "", builder->generateAdditionalResults(), req.id,
        BlackCalibrationHelper::RelativePriceError, allowChangingFallbacks, allowModelFallbacks, dontCalibrate);
    builder->engineFactory()->modelBuilders().insert(std::make_pair(req.id, fmmBuilder));
    return fmmBuilder;
}

QuantExt::FmmLsmEngineConfig fmmLsmEngineConfig(const EngineBuilder* builder,
                                                const QuantLib::ext::shared_ptr<FmmBuilder>& fmmBuilder,
                                                const bool observeIrModel) {
    FmmLsmEngineConfig cfg;
    cfg.lsm.trainingPaths = static_cast<Size>(parseInteger(builder->engineParameter("TrainingPaths", {}, false, "16384")));
    cfg.lsm.valuationPaths =
        static_cast<Size>(parseInteger(builder->engineParameter("ValuationPaths", {}, false, "16384")));
    cfg.lsm.trainingSeed = static_cast<BigNatural>(parseInteger(builder->engineParameter("TrainingSeed", {}, false, "42")));
    cfg.lsm.valuationSeed =
        static_cast<BigNatural>(parseInteger(builder->engineParameter("ValuationSeed", {}, false, "424242")));
    cfg.lsm.trainingSequence =
        parseSequenceType(builder->engineParameter("TrainingSequence", {}, false, "MersenneTwister"));
    cfg.lsm.valuationSequence =
        parseSequenceType(builder->engineParameter("ValuationSequence", {}, false, "MersenneTwister"));
    cfg.lsm.basisOrder = static_cast<Size>(parseInteger(builder->engineParameter("BasisOrder", {}, false, "2")));
    cfg.dualBound = parseBool(builder->engineParameter("DualBound", {}, false, "false"));
    cfg.dualOuterPaths = static_cast<Size>(parseInteger(builder->engineParameter("DualOuterPaths", {}, false, "512")));
    cfg.dualInnerPaths = static_cast<Size>(parseInteger(builder->engineParameter("DualInnerPaths", {}, false, "64")));
    cfg.dualSeed = static_cast<BigNatural>(parseInteger(builder->engineParameter("DualSeed", {}, false, "20260920")));
    const std::string policyMode = builder->engineParameter("PolicyMode", {}, false, "Retrain");
    QL_REQUIRE(policyMode == "Retrain" || policyMode == "Frozen",
               "FMM engine parameter PolicyMode must be Retrain or Frozen, got " << policyMode);
    cfg.policyMode = policyMode == "Frozen" ? FmmPolicyMode::Frozen : FmmPolicyMode::Retrain;
    cfg.gridToleranceDays =
        static_cast<Natural>(parseInteger(builder->modelParameter("GridToleranceDays", {}, false, "3")));
    // the calibrated IrModel adapter notifies after each recalibration; the cross-currency engine
    // calibrates lazily at pricing time and observes the builder itself instead
    if (observeIrModel)
        cfg.observables = {fmmBuilder->irModel()};
    else
        cfg.observables = {fmmBuilder};
    if (builder->generateAdditionalResults())
        cfg.calibrationResults = [fmmBuilder]() { return fmmBuilder->calibrationInfo().additionalResults(); };
    return cfg;
}

QuantLib::ext::shared_ptr<PricingEngine> FmmLsmSwaptionEngineBuilder::engineImpl(
    const string& id, const std::vector<string>& keys, const std::vector<Date>& dates,
    const std::vector<Date>& maturities, const std::vector<std::vector<Real>>& strikes,
    const std::vector<std::vector<Real>>& fxStrikes, const bool isAmerican, const std::string& discountCurve,
    const std::string& securitySpread, const CamOrLgmModel& modelOverwrite) {
    DLOG("Building FMM LSM Swaption engine for trade " << id);
    QL_REQUIRE(!isAmerican, "FmmLsmSwaptionEngineBuilder: American exercise is not supported");
    QL_REQUIRE(std::holds_alternative<std::monostate>(modelOverwrite),
               "FmmLsmSwaptionEngineBuilder: an externally supplied CAM / LGM model cannot be used with the FMM");
    QL_REQUIRE(discountCurve.empty() && securitySpread.empty(),
               "FmmLsmSwaptionEngineBuilder: a separate discount curve / security spread is not supported "
               "(single-curve baseline: the swap index discount curve is the model curve)");
    if (keys.size() == 2)
        return fmmXccyReducedEngine(this, id, keys, dates, maturities, strikes);
    QL_REQUIRE(keys.size() == 1, "FmmLsmSwaptionEngineBuilder: " << keys.size() << " currency keys are not supported");
    FmmModelRequest req;
    req.id = id;
    req.key = keys.front();
    req.expiries = dates;
    req.maturities = maturities;
    req.strikes = strikes.front();
    std::set<Date> gridDates(dates.begin(), dates.end());
    gridDates.insert(maturities.begin(), maturities.end());
    req.gridDates.assign(gridDates.begin(), gridDates.end());
    auto fmmBuilder = fmmModelBuilder(this, req);
    return QuantLib::ext::make_shared<FmmLsmMultiLegOptionEngine>(fmmBuilder->modelAsFmm(), fmmBuilder->grid(),
                                                                   fmmLsmEngineConfig(this, fmmBuilder));
}

QuantLib::ext::shared_ptr<QuantExt::PricingEngine> FmmLsmCallableBondEngineBuilder::engineImpl(
    const std::string& id, const std::string& ccy, const std::string& creditCurveId, const std::string& securityId,
    const std::string& referenceCurveId, const std::string& incomeCurveId, const QuantLib::Date& maturityDate) {
    DLOG("Building FMM LSM CallableBond engine for trade " << id);
    const std::string marketConfig = configuration(MarketContext::pricing);
    QL_REQUIRE(incomeCurveId.empty() || incomeCurveId == referenceCurveId,
               "FmmLsmCallableBondEngineBuilder: a separate income curve is not supported in this release");
    // coterminal basket on the reference calibration grid dates before maturity (as the LGM
    // callable bond builders); the model grid holds those dates plus the maturity
    const std::string referenceCalibrationGrid = modelParameter("ReferenceCalibrationGrid", {}, true, "");
    // unadjusted anniversaries (NullCalendar): the model grid must carry the bond's coupon / call
    // anniversaries, which the builder interface does not expose; the default TARGET adjustment
    // would move a weekend anniversary and shift the mapped pay dates by a day
    DateGrid grid(referenceCalibrationGrid, NullCalendar());
    FmmModelRequest req;
    req.id = id;
    req.key = ccy;
    for (const Date& d : grid.dates())
        if (d < maturityDate) {
            req.expiries.push_back(d);
            req.maturities.push_back(maturityDate);
            req.strikes.push_back(Null<Real>());
        }
    std::set<Date> gridDates(req.expiries.begin(), req.expiries.end());
    gridDates.insert(maturityDate);
    req.gridDates.assign(gridDates.begin(), gridDates.end());
    auto fmmBuilder = fmmModelBuilder(this, req);

    Handle<YieldTermStructure> referenceCurve = market_->yieldCurve(referenceCurveId, marketConfig);
    Handle<Quote> spread;
    if (!securityId.empty()) {
        try {
            spread = market_->securitySpread(securityId, marketConfig);
        } catch (const std::exception& e) {
            DLOG("FmmLsmCallableBondEngineBuilder: no security spread for " << securityId << " (" << e.what()
                                                                              << "), pricing without spread");
        }
    }
    Handle<DefaultProbabilityTermStructure> creditCurve;
    if (!creditCurveId.empty())
        creditCurve = market_->defaultCurve(creditCurveId, marketConfig)->curve(); // rejected by the engine
    return QuantLib::ext::make_shared<FmmLsmCallableBondEngine>(fmmBuilder->modelAsFmm(), fmmBuilder->grid(),
                                                                 fmmLsmEngineConfig(this, fmmBuilder), referenceCurve,
                                                                 spread, creditCurve);
}

} // namespace data
} // namespace ore
