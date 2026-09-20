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

#include <ored/model/fmmdata.hpp>
#include <ored/utilities/log.hpp>
#include <ored/utilities/parsers.hpp>
#include <ored/utilities/to_string.hpp>

#include <boost/algorithm/string.hpp>

namespace ore {
namespace data {

QuantExt::FmmParametrization::LocalVolType parseFmmVolatilityType(const std::string& s) {
    const std::string u = boost::algorithm::to_upper_copy(s);
    if (u == "DISPLACEDDIFFUSION" || u == "DD" || u == "SHIFTEDLOGNORMAL")
        return QuantExt::FmmParametrization::LocalVolType::DisplacedDiffusion;
    if (u == "NORMAL")
        return QuantExt::FmmParametrization::LocalVolType::Normal;
    QL_FAIL("FMM VolatilityType '" << s << "' not recognized (DisplacedDiffusion, Normal)");
}

QuantExt::FmmSwaptionApproxMethod parseFmmApproximationMethod(const std::string& s) {
    const std::string u = boost::algorithm::to_upper_copy(s);
    if (u == "STRIKEFLATNORMAL" || u == "STRIKEFLAT")
        return QuantExt::FmmSwaptionApproxMethod::StrikeFlatNormal;
    if (u == "EFFECTIVESHIFT")
        return QuantExt::FmmSwaptionApproxMethod::EffectiveShift;
    QL_FAIL("FMM ApproximationMethod '" << s << "' not recognized (StrikeFlatNormal, EffectiveShift)");
}

Real parseFmmShift(const std::string& s) {
    const std::string u = boost::algorithm::to_upper_copy(boost::algorithm::trim_copy(s));
    if (u == "1/TAU" || u == "HJM" || u == "REPLICATION")
        return Null<Real>();
    return parseReal(s);
}

std::string to_string(const QuantExt::FmmParametrization::LocalVolType& t) {
    return t == QuantExt::FmmParametrization::LocalVolType::Normal ? "Normal" : "DisplacedDiffusion";
}

std::string to_string(const QuantExt::FmmSwaptionApproxMethod& m) {
    return m == QuantExt::FmmSwaptionApproxMethod::EffectiveShift ? "EffectiveShift" : "StrikeFlatNormal";
}

FmmData::FmmData() : IrModelData("FMM", "", CalibrationType::None) { reset(); }

void FmmData::clear() {
    optionExpiries_.clear();
    optionTerms_.clear();
    optionStrikes_.clear();
    gridDates_.clear();
}

void FmmData::reset() {
    IrModelData::reset();
    grid_ = 3 * Months;
    shift_ = Null<Real>();
    volType_ = QuantExt::FmmParametrization::LocalVolType::DisplacedDiffusion;
    rhoInf_ = 0.6;
    beta_ = 0.08;
    factors_ = 3;
    calibrateVol_ = true;
    volParamType_ = ParamType::Piecewise;
    volTimes_ = {};
    volValues_ = {0.0025};
    approx_ = QuantExt::FmmSwaptionApproxMethod::EffectiveShift;
    mc_ = McCorrection();
    subSteps_ = 1;
    gridToleranceDays_ = 3;
}

void FmmData::validate() const {
    QL_REQUIRE(grid_.length() > 0, "FmmData: Grid tenor must be positive");
    QL_REQUIRE(shift_ == Null<Real>() || shift_ >= 0.0, "FmmData: Shift must be non-negative or 1/tau");
    QL_REQUIRE(rhoInf_ > -1.0 && rhoInf_ <= 1.0, "FmmData: RhoInf must be in (-1, 1]");
    QL_REQUIRE(beta_ >= 0.0, "FmmData: Beta must be non-negative");
    QL_REQUIRE(factors_ >= 1, "FmmData: Factors must be at least 1");
    QL_REQUIRE(!volValues_.empty(), "FmmData: Volatility needs at least one initial value");
    if (volParamType_ == ParamType::Constant)
        QL_REQUIRE(volTimes_.empty() && volValues_.size() == 1,
                   "FmmData: constant volatility needs an empty time grid and one initial value");
    else
        QL_REQUIRE(volValues_.size() == 1 || volValues_.size() == volTimes_.size() + 1,
                   "FmmData: piecewise volatility needs one initial value or one per segment ("
                       << volTimes_.size() + 1 << ")");
    for (const Real v : volValues_)
        QL_REQUIRE(v > 0.0, "FmmData: volatility initial values must be positive");
    QL_REQUIRE(subSteps_ >= 1, "FmmData: SubSteps must be at least 1");
    if (mc_.enabled) {
        QL_REQUIRE(mc_.reps >= 2 && mc_.pathsPerRep > 0, "FmmData: McCorrected needs Reps >= 2 and PathsPerRep > 0");
        QL_REQUIRE(mc_.toleranceBp > 0.0 && mc_.maxIterations >= 1, "FmmData: McCorrected tolerance/iterations invalid");
    }
    QL_REQUIRE(optionExpiries_.size() == optionTerms_.size() && optionStrikes_.size() == optionExpiries_.size(),
               "FmmData: calibration swaption expiries/terms/strikes sizes differ");
    if (calibrationType_ != CalibrationType::None)
        QL_REQUIRE(calibrationType_ == CalibrationType::Bootstrap || calibrationType_ == CalibrationType::BestFit,
                   "FmmData: CalibrationType must be None, Bootstrap or BestFit");
    if (calibrateVol_ && calibrationType_ != CalibrationType::None)
        QL_REQUIRE(!optionExpiries_.empty(), "FmmData: volatility calibration needs calibration swaptions");
}

void FmmData::fromXML(XMLNode* node) {
    reset();
    qualifier_ = XMLUtils::getAttribute(node, "key");
    LOG("FMM with attribute (key) = " << qualifier_);

    IrModelData::fromXML(node); // CalibrationType

    grid_ = parsePeriod(XMLUtils::getChildValue(node, "Grid", false, "3M"));
    shift_ = parseFmmShift(XMLUtils::getChildValue(node, "Shift", false, "1/tau"));
    volType_ = parseFmmVolatilityType(XMLUtils::getChildValue(node, "VolatilityType", false, "DisplacedDiffusion"));
    rhoInf_ = XMLUtils::getChildValueAsDouble(node, "RhoInf", false, 0.6);
    beta_ = XMLUtils::getChildValueAsDouble(node, "Beta", false, 0.08);
    factors_ = static_cast<Size>(XMLUtils::getChildValueAsInt(node, "Factors", false, 3));

    if (XMLNode* volNode = XMLUtils::getChildNode(node, "Volatility")) {
        calibrateVol_ = XMLUtils::getChildValueAsBool(volNode, "Calibrate", false, true);
        volParamType_ = parseParamType(XMLUtils::getChildValue(volNode, "ParamType", false, "Piecewise"));
        volTimes_ = XMLUtils::getChildrenValuesAsDoublesCompact(volNode, "TimeGrid", false);
        volValues_ = XMLUtils::getChildrenValuesAsDoublesCompact(volNode, "InitialValue", true);
    }

    approx_ = parseFmmApproximationMethod(XMLUtils::getChildValue(node, "ApproximationMethod", false, "EffectiveShift"));

    if (XMLNode* mcNode = XMLUtils::getChildNode(node, "McCorrected")) {
        mc_.enabled = XMLUtils::getChildValueAsBool(mcNode, "Enabled", false, true);
        mc_.pathsPerRep = static_cast<Size>(XMLUtils::getChildValueAsInt(mcNode, "PathsPerRep", false, 16384));
        mc_.reps = static_cast<Size>(XMLUtils::getChildValueAsInt(mcNode, "Reps", false, 8));
        mc_.seed = static_cast<BigNatural>(XMLUtils::getChildValueAsInt(mcNode, "Seed", false, 9100));
        mc_.maxIterations = static_cast<Size>(XMLUtils::getChildValueAsInt(mcNode, "MaxIterations", false, 3));
        mc_.toleranceBp = XMLUtils::getChildValueAsDouble(mcNode, "ToleranceBp", false, 0.2);
    }

    subSteps_ = static_cast<Size>(XMLUtils::getChildValueAsInt(node, "SubSteps", false, 1));
    gridToleranceDays_ = static_cast<Natural>(XMLUtils::getChildValueAsInt(node, "GridToleranceDays", false, 3));

    if (XMLNode* optionsNode = XMLUtils::getChildNode(node, "CalibrationSwaptions")) {
        optionExpiries() = XMLUtils::getChildrenValuesAsStrings(optionsNode, "Expiries", false);
        optionTerms() = XMLUtils::getChildrenValuesAsStrings(optionsNode, "Terms", false);
        QL_REQUIRE(optionExpiries().size() == optionTerms().size(),
                   "FmmData: vector size mismatch in swaption expiries/terms for key " << qualifier_);
        optionStrikes() = XMLUtils::getChildrenValuesAsStrings(optionsNode, "Strikes", false);
        if (!optionStrikes().empty())
            QL_REQUIRE(optionStrikes().size() == optionExpiries().size(),
                       "FmmData: vector size mismatch in swaption expiries/strikes for key " << qualifier_);
        else
            optionStrikes().resize(optionExpiries().size(), "ATM");
        for (Size i = 0; i < optionExpiries().size(); i++)
            LOG("FMM calibration swaption " << optionExpiries()[i] << " x " << optionTerms()[i] << " "
                                            << optionStrikes()[i]);
    }
    validate();
    LOG("FmmData done");
}

XMLNode* FmmData::toXML(XMLDocument& doc) const {
    XMLNode* node = IrModelData::toXML(doc);
    XMLUtils::addAttribute(doc, node, "key", qualifier_);
    XMLUtils::addChild(doc, node, "Grid", ore::data::to_string(grid_));
    XMLUtils::addChild(doc, node, "Shift", shift_ == Null<Real>() ? std::string("1/tau") : ore::data::to_string(shift_));
    XMLUtils::addChild(doc, node, "VolatilityType", to_string(volType_));
    XMLUtils::addChild(doc, node, "RhoInf", rhoInf_);
    XMLUtils::addChild(doc, node, "Beta", beta_);
    XMLUtils::addChild(doc, node, "Factors", static_cast<int>(factors_));
    XMLNode* volNode = XMLUtils::addChild(doc, node, "Volatility");
    XMLUtils::addChild(doc, volNode, "Calibrate", calibrateVol_);
    XMLUtils::addGenericChild(doc, volNode, "ParamType", volParamType_);
    XMLUtils::addGenericChildAsList(doc, volNode, "TimeGrid", volTimes_);
    XMLUtils::addGenericChildAsList(doc, volNode, "InitialValue", volValues_);
    XMLUtils::addChild(doc, node, "ApproximationMethod", to_string(approx_));
    XMLNode* mcNode = XMLUtils::addChild(doc, node, "McCorrected");
    XMLUtils::addChild(doc, mcNode, "Enabled", mc_.enabled);
    XMLUtils::addChild(doc, mcNode, "PathsPerRep", static_cast<int>(mc_.pathsPerRep));
    XMLUtils::addChild(doc, mcNode, "Reps", static_cast<int>(mc_.reps));
    XMLUtils::addChild(doc, mcNode, "Seed", static_cast<int>(mc_.seed));
    XMLUtils::addChild(doc, mcNode, "MaxIterations", static_cast<int>(mc_.maxIterations));
    XMLUtils::addChild(doc, mcNode, "ToleranceBp", mc_.toleranceBp);
    XMLUtils::addChild(doc, node, "SubSteps", static_cast<int>(subSteps_));
    XMLUtils::addChild(doc, node, "GridToleranceDays", static_cast<int>(gridToleranceDays_));
    XMLNode* calibrationSwaptionsNode = XMLUtils::addChild(doc, node, "CalibrationSwaptions");
    XMLUtils::addGenericChildAsList(doc, calibrationSwaptionsNode, "Expiries", optionExpiries());
    XMLUtils::addGenericChildAsList(doc, calibrationSwaptionsNode, "Terms", optionTerms());
    XMLUtils::addGenericChildAsList(doc, calibrationSwaptionsNode, "Strikes", optionStrikes());
    return node;
}

} // namespace data
} // namespace ore
