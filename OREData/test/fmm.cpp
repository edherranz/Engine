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

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <sstream>

#include <ored/configuration/conventions.hpp>
#include <ored/marketdata/marketimpl.hpp>
#include <ored/model/fmmbuilder.hpp>
#include <ored/model/fmmdata.hpp>
#include <ored/portfolio/builders/bond.hpp>
#include <ored/portfolio/builders/callablebond.hpp>
#include <ored/portfolio/builders/fmm.hpp>
#include <ored/portfolio/builders/swap.hpp>
#include <ored/portfolio/builders/swaption.hpp>
#include <ored/portfolio/callablebond.hpp>
#include <ored/portfolio/enginedata.hpp>
#include <ored/portfolio/swaption.hpp>
#include <ored/utilities/indexparser.hpp>
#include <ored/utilities/to_string.hpp>
#include <ored/utilities/toplevelfixture.hpp>

#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/volatility/optionlet/constantoptionletvol.hpp>
#include <ql/termstructures/volatility/swaption/swaptionconstantvol.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

using namespace QuantLib;
using namespace QuantExt;
using namespace boost::unit_test_framework;
using namespace std;
using namespace ore::data;

namespace {

// flat USD-SOFR market with an OIS swap index family, a normal swaption vol and a security spread,
// all behind quotes so that the recalibration modes can be exercised
class FmmTestMarket : public MarketImpl {
public:
    FmmTestMarket(const Date& asof, const Real flatRate, const Real normalVol) : MarketImpl(false) {
        asof_ = asof;
        curveQuote_ = QuantLib::ext::make_shared<SimpleQuote>(flatRate);
        volQuote_ = QuantLib::ext::make_shared<SimpleQuote>(normalVol);
        spreadQuote_ = QuantLib::ext::make_shared<SimpleQuote>(0.0);
        Handle<YieldTermStructure> yts(QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), Handle<Quote>(curveQuote_),
                                                                              Actual365Fixed()));
        Handle<IborIndex> sofr(parseIborIndex("USD-SOFR", yts));
        iborIndices_[make_pair(Market::defaultConfiguration, "USD-SOFR")] = sofr;
        // recent fixings: an overnight period starting on a weekend uses the preceding business day's fixing
        for (Date d = asof - 10; d < asof; ++d)
            if (sofr->isValidFixingDate(d))
                sofr->addFixing(d, flatRate);

        auto conventions = QuantLib::ext::make_shared<Conventions>();
        conventions->add(QuantLib::ext::make_shared<OisConvention>("USD-SOFR-OIS", "0", "USD-SOFR", "A365", "NullCalendar",
                                                                   "0", "false", "Annual", "F", "F", "Forward"));
        conventions->add(QuantLib::ext::make_shared<SwapIndexConvention>("USD-CMS-2Y", "USD-SOFR-OIS"));
        conventions->add(QuantLib::ext::make_shared<SwapIndexConvention>("USD-CMS-30Y", "USD-SOFR-OIS"));
        InstrumentConventions::instance().setConventions(conventions);

        yieldCurves_[make_tuple(Market::defaultConfiguration, YieldCurveType::Discount, "USD")] = yts;
        yieldCurves_[make_tuple(Market::defaultConfiguration, YieldCurveType::Yield, "USD-SOFR")] = yts;
        addSwapIndex("USD-CMS-2Y", "USD-SOFR", Market::defaultConfiguration);
        addSwapIndex("USD-CMS-30Y", "USD-SOFR", Market::defaultConfiguration);
        swaptionIndexBases_[make_pair(Market::defaultConfiguration, "USD-SOFR")] = make_pair("USD-CMS-2Y", "USD-CMS-30Y");
        swaptionIndexBases_[make_pair(Market::defaultConfiguration, "USD")] = make_pair("USD-CMS-2Y", "USD-CMS-30Y");
        Handle<QuantLib::SwaptionVolatilityStructure> svs(QuantLib::ext::make_shared<QuantLib::ConstantSwaptionVolatility>(
            0, NullCalendar(), Following, Handle<Quote>(volQuote_), Actual365Fixed(), Normal, 0.0));
        swaptionCurves_[make_pair(Market::defaultConfiguration, "USD")] = svs;
        swaptionCurves_[make_pair(Market::defaultConfiguration, "USD-SOFR")] = svs;
        securitySpreads_[make_pair(Market::defaultConfiguration, "SEC-FMM")] = Handle<Quote>(spreadQuote_);
        Handle<QuantLib::OptionletVolatilityStructure> ovs(QuantLib::ext::make_shared<QuantLib::ConstantOptionletVolatility>(
            0, NullCalendar(), Following, Handle<Quote>(volQuote_), Actual365Fixed(), Normal, 0.0));
        capFloorCurves_[make_pair(Market::defaultConfiguration, "USD")] = ovs;
        capFloorCurves_[make_pair(Market::defaultConfiguration, "USD-SOFR")] = ovs;
    }
    QuantLib::ext::shared_ptr<SimpleQuote> curveQuote_, volQuote_, spreadQuote_;
};

const Date kAsof(19, September, 2026);

QuantLib::ext::shared_ptr<EngineData> fmmEngineData(const std::string& product, const bool dualBound,
                                                    const std::map<std::string, std::string>& extraModelParams = {}) {
    auto ed = QuantLib::ext::make_shared<EngineData>();
    ed->model(product) = "FMM";
    ed->engine(product) = "LSM";
    ed->modelParameters(product) = {{"Calibration", "Bootstrap"},
                                    {"CalibrationStrategy", "CoterminalATM"},
                                    {"Tolerance", "0.01"},
                                    {"Volatility", "0.0025"},
                                    {"VolatilityType", "DisplacedDiffusion"},
                                    {"Shift", "1/tau"},
                                    {"RhoInf", "0.6"},
                                    {"Beta", "0.08"},
                                    {"Factors", "3"},
                                    {"ApproximationMethod", "EffectiveShift"},
                                    {"Grid", "3M"},
                                    {"GridToleranceDays", "3"},
                                    {"ReferenceCalibrationGrid", "1Y,2Y,3Y,4Y"}};
    for (const auto& kv : extraModelParams)
        ed->modelParameters(product)[kv.first] = kv.second;
    ed->engineParameters(product) = {{"TrainingPaths", "16384"}, {"ValuationPaths", "16384"},
                                     {"TrainingSeed", "42"},     {"ValuationSeed", "424242"},
                                     {"BasisOrder", "2"},        {"DualBound", dualBound ? "true" : "false"},
                                     {"DualOuterPaths", "256"},  {"DualInnerPaths", "32"}};
    ed->globalParameters()["GenerateAdditionalResults"] = "true";
    ed->model("Swap") = "DiscountedCashflows";
    ed->engine("Swap") = "DiscountingSwapEngine";
    ed->model("Bond") = "DiscountedCashflows";
    ed->engine("Bond") = "DiscountingRiskyBondEngine";
    ed->engineParameters("Bond") = {{"TimestepPeriod", "6M"}};
    return ed;
}

QuantLib::ext::shared_ptr<EngineData> lgmEngineData(const std::string& product) {
    auto ed = QuantLib::ext::make_shared<EngineData>();
    ed->model(product) = "LGM";
    ed->engine(product) = "Grid";
    ed->modelParameters(product) = {{"Calibration", "Bootstrap"},
                                    {"CalibrationStrategy", "CoterminalATM"},
                                    {"Reversion", "0.02"},
                                    {"Volatility", "0.01"},
                                    {"VolatilityType", "Hagan"},
                                    {"ReversionType", "HullWhite"},
                                    {"Tolerance", "0.0001"},
                                    {"ShiftHorizon", "0.5"},
                                    {"ReferenceCalibrationGrid", "1Y,2Y,3Y,4Y"}};
    ed->engineParameters(product) = {{"sy", "3.0"}, {"ny", "10"}, {"sx", "3.0"}, {"nx", "10"}};
    ed->globalParameters()["GenerateAdditionalResults"] = "true";
    ed->model("Swap") = "DiscountedCashflows";
    ed->engine("Swap") = "DiscountingSwapEngine";
    ed->model("Bond") = "DiscountedCashflows";
    ed->engine("Bond") = "DiscountingRiskyBondEngine";
    ed->engineParameters("Bond") = {{"TimestepPeriod", "6M"}};
    return ed;
}

QuantLib::ext::shared_ptr<ore::data::Swaption>
bermudanPayerSwaption(const Real fixedRate, const std::string& noticePeriod = "", const std::string& start = "2026-09-19",
                      const std::string& end = "2031-09-19",
                      const std::vector<std::string>& exerciseDates = {"2027-09-19", "2028-09-19", "2029-09-19",
                                                                       "2030-09-19"}) {
    ScheduleData fixedSchedule(ScheduleRules(start, end, "1Y", "NullCalendar", "F", "F", "Forward"));
    ScheduleData floatSchedule(ScheduleRules(start, end, "3M", "NullCalendar", "F", "F", "Forward"));
    LegData fixedLeg(QuantLib::ext::make_shared<FixedLegData>(std::vector<Real>(1, fixedRate)), true, "USD", fixedSchedule,
                     "A365", std::vector<Real>(1, 1.0));
    LegData floatLeg(QuantLib::ext::make_shared<FloatingLegData>("USD-SOFR", 0, true, std::vector<Real>(1, 0.0)), false,
                     "USD", floatSchedule, "A365", std::vector<Real>(1, 1.0));
    Envelope env("CP1");
    OptionData optionData("Long", "Call", "Bermudan", true, exerciseDates, "Physical", "", PremiumData(), {}, {},
                          noticePeriod, noticePeriod.empty() ? "" : "NullCalendar",
                          noticePeriod.empty() ? "" : "Unadjusted");
    return QuantLib::ext::make_shared<ore::data::Swaption>(env, optionData, std::vector<LegData>{fixedLeg, floatLeg});
}

QuantLib::ext::shared_ptr<ore::data::CallableBond> callableBondFromXml(const std::string& extraCallData = "") {
    const std::string xml =
        "<Trade id=\"CB-FMM\"><TradeType>CallableBond</TradeType>"
        "<Envelope><CounterParty>CP1</CounterParty><NettingSetId>NS1</NettingSetId></Envelope>"
        "<CallableBondData><BondData>"
        "<IssuerId>ISSUER</IssuerId><CreditCurveId></CreditCurveId><CreditRisk>false</CreditRisk>"
        "<SecurityId>SEC-FMM</SecurityId><ReferenceCurveId>USD-SOFR</ReferenceCurveId>"
        "<SettlementDays>0</SettlementDays><Calendar>NullCalendar</Calendar><IssueDate>2026-09-19</IssueDate>"
        "<LegData><LegType>Fixed</LegType><Payer>false</Payer><Currency>USD</Currency>"
        "<Notionals><Notional>1</Notional></Notionals><DayCounter>A365</DayCounter><PaymentConvention>F</PaymentConvention>"
        "<FixedLegData><Rates><Rate>0.031</Rate></Rates></FixedLegData>"
        "<ScheduleData><Rules><StartDate>2026-09-19</StartDate><EndDate>2031-09-19</EndDate><Tenor>1Y</Tenor>"
        "<Calendar>NullCalendar</Calendar><Convention>F</Convention><TermConvention>F</TermConvention><Rule>Forward</Rule>"
        "</Rules></ScheduleData></LegData></BondData>"
        "<CallData><Styles><Style>Bermudan</Style></Styles>"
        "<ScheduleData><Dates><Dates><Date>2027-09-19</Date><Date>2028-09-19</Date><Date>2029-09-19</Date>"
        "<Date>2030-09-19</Date></Dates></Dates></ScheduleData>"
        "<Prices><Price>1.0</Price></Prices><PriceTypes><PriceType>Dirty</PriceType></PriceTypes>"
        "<IncludeAccruals><IncludeAccrual>true</IncludeAccrual></IncludeAccruals>" +
        extraCallData + "</CallData></CallableBondData></Trade>";
    XMLDocument doc;
    doc.fromXMLString(xml);
    auto trade = QuantLib::ext::make_shared<ore::data::CallableBond>();
    trade->fromXML(doc.getFirstNode("Trade"));
    return trade;
}

} // namespace

BOOST_FIXTURE_TEST_SUITE(OREDataTestSuite, ore::data::TopLevelFixture)

BOOST_AUTO_TEST_SUITE(FmmTests)

BOOST_AUTO_TEST_CASE(testFmmDataXmlRoundTrip) {
    BOOST_TEST_MESSAGE("Testing FmmData XML parsing, serialisation round trip and validation...");
    const std::string xml =
        "<FMM key=\"USD-SOFR\"><CalibrationType>Bootstrap</CalibrationType><Grid>3M</Grid><Shift>1/tau</Shift>"
        "<VolatilityType>DisplacedDiffusion</VolatilityType><RhoInf>0.7</RhoInf><Beta>0.05</Beta><Factors>2</Factors>"
        "<Volatility><Calibrate>Y</Calibrate><ParamType>Piecewise</ParamType><TimeGrid>1.0,2.0</TimeGrid>"
        "<InitialValue>0.003</InitialValue></Volatility><ApproximationMethod>EffectiveShift</ApproximationMethod>"
        "<McCorrected><Enabled>true</Enabled><PathsPerRep>8192</PathsPerRep><Reps>4</Reps><Seed>7</Seed>"
        "<MaxIterations>2</MaxIterations><ToleranceBp>0.5</ToleranceBp></McCorrected><SubSteps>2</SubSteps>"
        "<GridToleranceDays>4</GridToleranceDays><CalibrationSwaptions><Expiries>1Y,2Y,3Y</Expiries>"
        "<Terms>4Y,3Y,2Y</Terms><Strikes>ATM,ATM,0.03</Strikes></CalibrationSwaptions></FMM>";
    XMLDocument doc;
    doc.fromXMLString(xml);
    FmmData data;
    data.fromXML(doc.getFirstNode("FMM"));
    BOOST_CHECK_EQUAL(data.qualifier(), "USD-SOFR");
    BOOST_CHECK(data.calibrationType() == CalibrationType::Bootstrap);
    BOOST_CHECK(data.grid() == 3 * Months);
    BOOST_CHECK(data.shift() == Null<Real>());
    BOOST_CHECK(data.volatilityType() == FmmParametrization::LocalVolType::DisplacedDiffusion);
    BOOST_CHECK_CLOSE(data.rhoInf(), 0.7, 1e-12);
    BOOST_CHECK_CLOSE(data.beta(), 0.05, 1e-12);
    BOOST_CHECK_EQUAL(data.factors(), Size(2));
    BOOST_CHECK(data.calibrateVolatility());
    BOOST_CHECK_EQUAL(data.volTimes().size(), Size(2));
    BOOST_CHECK_EQUAL(data.volValues().size(), Size(1));
    BOOST_CHECK(data.approximationMethod() == FmmSwaptionApproxMethod::EffectiveShift);
    BOOST_CHECK(data.mcCorrection().enabled);
    BOOST_CHECK_EQUAL(data.mcCorrection().reps, Size(4));
    BOOST_CHECK_EQUAL(data.subSteps(), Size(2));
    BOOST_CHECK_EQUAL(data.gridToleranceDays(), Natural(4));
    BOOST_CHECK_EQUAL(data.optionExpiries().size(), Size(3));
    BOOST_CHECK_EQUAL(data.optionStrikes()[2], "0.03");
    // round trip
    XMLDocument out;
    XMLNode* node = data.toXML(out);
    FmmData again;
    again.fromXML(node);
    BOOST_CHECK_EQUAL(again.qualifier(), data.qualifier());
    BOOST_CHECK(again.grid() == data.grid());
    BOOST_CHECK_CLOSE(again.rhoInf(), data.rhoInf(), 1e-12);
    BOOST_CHECK_EQUAL(again.factors(), data.factors());
    BOOST_CHECK_EQUAL(again.volTimes().size(), data.volTimes().size());
    BOOST_CHECK(again.mcCorrection().enabled == data.mcCorrection().enabled);
    BOOST_CHECK_EQUAL(again.mcCorrection().pathsPerRep, data.mcCorrection().pathsPerRep);
    BOOST_CHECK_EQUAL(again.optionTerms().size(), Size(3));
    // invalid combinations are rejected
    FmmData bad = data;
    bad.factors() = 0;
    BOOST_CHECK_THROW(bad.validate(), QuantLib::Error);
    FmmData bad2 = data;
    bad2.volValues() = {-0.1};
    BOOST_CHECK_THROW(bad2.validate(), QuantLib::Error);
    FmmData bad3 = data;
    bad3.calibrationType() = CalibrationType::FirstBestFitThenBootstrap;
    BOOST_CHECK_THROW(bad3.validate(), QuantLib::Error);
}

BOOST_AUTO_TEST_CASE(testFmmBermudanSwaptionThroughEngineFactory) {
    BOOST_TEST_MESSAGE("Testing a Bermudan swaption priced through the engine factory with model FMM / engine LSM "
                       "(coterminal ATM bootstrap) against ORE's LGM grid engine on the same trade...");
    Settings::instance().evaluationDate() = kAsof;
    auto market = QuantLib::ext::make_shared<FmmTestMarket>(kAsof, 0.03, 0.0080);

    auto fmmFactory = QuantLib::ext::make_shared<EngineFactory>(fmmEngineData("BermudanSwaption", true), market);
    auto trade = bermudanPayerSwaption(0.031);
    trade->build(fmmFactory);
    const Real npvFmm = trade->instrument()->NPV();
    const auto& add = trade->instrument()->additionalResults();
    BOOST_TEST_MESSAGE("FMM/LSM Bermudan payer NPV " << npvFmm);
    BOOST_REQUIRE(add.count("fmmLsmLowerBound") > 0);
    const Real se = QuantLib::ext::any_cast<Real>(add.at("fmmLsmLowerBoundStdError"));
    BOOST_CHECK_CLOSE(QuantLib::ext::any_cast<Real>(add.at("fmmLsmLowerBound")), npvFmm, 1e-10);
    BOOST_REQUIRE(add.count("fmmDualityGap") > 0);
    BOOST_CHECK(QuantLib::ext::any_cast<Real>(add.at("fmmDualityGap")) >
                -3.0 * QuantLib::ext::any_cast<Real>(add.at("fmmDualityGapStdError")));
    BOOST_REQUIRE(add.count("calibrationError") > 0);
    BOOST_CHECK_EQUAL(QuantLib::ext::any_cast<std::string>(add.at("modelType")), "FMM");
    const auto modelVols = QuantLib::ext::any_cast<std::vector<Real>>(add.at("fmmCalibrationModelVols"));
    const auto marketVols = QuantLib::ext::any_cast<std::vector<Real>>(add.at("fmmCalibrationMarketVols"));
    BOOST_REQUIRE_EQUAL(modelVols.size(), Size(4));
    for (Size i = 0; i < modelVols.size(); ++i) {
        BOOST_TEST_MESSAGE("calibration swaption " << i << ": market " << marketVols[i] * 1e4 << " bp, model "
                                                   << modelVols[i] * 1e4 << " bp");
        BOOST_CHECK_SMALL(marketVols[i] - 0.0080, 1e-12);
        BOOST_CHECK_SMALL(modelVols[i] - marketVols[i], 1e-7); // analytic bootstrap (solver tolerance)
    }
    // the calibration helpers follow the SOFR fixing calendar while the trade uses NullCalendar dates
    // (weekend anniversaries): the grid mapping absorbs this within GridToleranceDays
    const Real mismatch = QuantLib::ext::any_cast<Real>(add.at("fmmGridMaxDateMismatchDays"));
    BOOST_TEST_MESSAGE("grid max date mismatch (basket vs trade dates): " << mismatch << " days");
    BOOST_CHECK(mismatch <= 3.0);
    BOOST_CHECK_EQUAL(QuantLib::ext::any_cast<std::string>(add.at("fmmExerciseProbabilityType")), "unconditional");
    BOOST_CHECK_EQUAL(fmmFactory->modelBuilders().size(), Size(1));

    auto lgmFactory = QuantLib::ext::make_shared<EngineFactory>(lgmEngineData("BermudanSwaption"), market);
    auto lgmTrade = bermudanPayerSwaption(0.031);
    lgmTrade->build(lgmFactory);
    const Real npvLgm = lgmTrade->instrument()->NPV();
    BOOST_TEST_MESSAGE("LGM/Grid Bermudan payer NPV " << npvLgm << "; FMM - LGM = " << npvFmm - npvLgm << " ("
                                                     << 100.0 * (npvFmm - npvLgm) / npvLgm << "% of value, LSM s.e. " << se
                                                     << ")");
    // different factor structures calibrated to the same coterminals: agreement at the few-percent
    // level is a model difference (reported), the plumbing check is the 10% band plus MC noise
    BOOST_CHECK_MESSAGE(std::fabs(npvFmm - npvLgm) < 0.10 * std::fabs(npvLgm) + 3.0 * se,
                        "FMM vs LGM Bermudan: " << npvFmm << " vs " << npvLgm);
}

BOOST_AUTO_TEST_CASE(testFmmRecalibrationModes) {
    BOOST_TEST_MESSAGE("Testing FmmBuilder lazy recalibration: vol bump + recalibrate, curve bump, and "
                       "newCalcWithoutRecalibration (frozen parameters)...");
    Settings::instance().evaluationDate() = kAsof;
    auto market = QuantLib::ext::make_shared<FmmTestMarket>(kAsof, 0.03, 0.0080);
    auto factory = QuantLib::ext::make_shared<EngineFactory>(fmmEngineData("BermudanSwaption", false), market);
    auto trade = bermudanPayerSwaption(0.031);
    trade->build(factory);
    const Real npv0 = trade->instrument()->NPV();
    const Real se = QuantLib::ext::any_cast<Real>(trade->instrument()->additionalResults().at("fmmLsmLowerBoundStdError"));
    BOOST_REQUIRE_EQUAL(factory->modelBuilders().size(), Size(1));
    auto builder = QuantLib::ext::dynamic_pointer_cast<FmmBuilder>(factory->modelBuilders().begin()->second);
    BOOST_REQUIRE(builder);
    BOOST_CHECK(!builder->requiresRecalibration());

    // vol bump: the builder notices, recalibration changes the value
    market->volQuote_->setValue(0.0100);
    BOOST_CHECK(builder->requiresRecalibration());
    builder->recalibrate();
    const Real npvVolUp = trade->instrument()->NPV();
    BOOST_TEST_MESSAGE("vol 80 -> 100 bp: NPV " << npv0 << " -> " << npvVolUp);
    BOOST_CHECK_MESSAGE(npvVolUp > npv0 + 3.0 * se, "higher vol must raise the Bermudan value beyond noise");
    const auto modelVols =
        QuantLib::ext::any_cast<std::vector<Real>>(trade->instrument()->additionalResults().at("fmmCalibrationModelVols"));
    for (const Real v : modelVols)
        BOOST_CHECK_SMALL(v - 0.0100, 1e-9);

    // frozen-parameter mode (scenario flow: the market moves, then newCalcWithoutRecalibration):
    // the vol move is ignored, calibrated parameters and value stay put
    market->volQuote_->setValue(0.0120);
    builder->newCalcWithoutRecalibration();
    BOOST_CHECK(!builder->requiresRecalibration());
    const Real npvFrozen = trade->instrument()->NPV();
    BOOST_TEST_MESSAGE("vol 100 -> 120 bp without recalibration: NPV " << npvFrozen << " (frozen at " << npvVolUp << ")");
    BOOST_CHECK_SMALL(npvFrozen - npvVolUp, 1e-12);
    const auto frozenVols =
        QuantLib::ext::any_cast<std::vector<Real>>(trade->instrument()->additionalResults().at("fmmCalibrationModelVols"));
    for (const Real v : frozenVols)
        BOOST_CHECK_SMALL(v - 0.0100, 1e-9);

    // next scenario: another vol move, this time followed by a recalibration
    market->volQuote_->setValue(0.0110);
    builder->recalibrate();
    const Real npvVolUp2 = trade->instrument()->NPV();
    const auto vols110 =
        QuantLib::ext::any_cast<std::vector<Real>>(trade->instrument()->additionalResults().at("fmmCalibrationModelVols"));
    for (const Real v : vols110)
        BOOST_CHECK_SMALL(v - 0.0110, 1e-9);
    BOOST_TEST_MESSAGE("vol 120 -> 110 bp with recalibration: NPV " << npvVolUp2);
    BOOST_CHECK_MESSAGE(npvVolUp2 > npvVolUp + 3.0 * se, "the 110 bp value must exceed the 100 bp value beyond noise");

    // curve bump: the engine observes the curve; recalibration keeps the vols, the value moves
    market->curveQuote_->setValue(0.035);
    builder->recalibrate();
    const Real npvCurveUp = trade->instrument()->NPV();
    BOOST_TEST_MESSAGE("curve 3% -> 3.5%: NPV " << npvVolUp2 << " -> " << npvCurveUp);
    BOOST_CHECK_MESSAGE(std::fabs(npvCurveUp - npvVolUp2) > 3.0 * se, "a 50 bp curve move must change the payer value");
}

BOOST_AUTO_TEST_CASE(testFmmCallableBondThroughEngineFactory) {
    BOOST_TEST_MESSAGE("Testing a callable bond priced through the engine factory with model FMM / engine LSM "
                       "against ORE's LGM grid callable bond engine, with and without an issuer spread...");
    Settings::instance().evaluationDate() = kAsof;
    auto market = QuantLib::ext::make_shared<FmmTestMarket>(kAsof, 0.03, 0.0080);

    auto fmmFactory = QuantLib::ext::make_shared<EngineFactory>(fmmEngineData("CallableBond", false), market);
    auto trade = callableBondFromXml();
    trade->build(fmmFactory);
    const Real npvFmm = trade->instrument()->NPV();
    const auto& add = trade->instrument()->additionalResults();
    BOOST_REQUIRE(add.count("callPutValue") > 0);
    const Real callFmm = QuantLib::ext::any_cast<Real>(add.at("callPutValue"));
    const Real stripped = QuantLib::ext::any_cast<Real>(add.at("strippedBondNpv"));
    const Real se = QuantLib::ext::any_cast<Real>(add.at("fmmLsmLowerBoundStdError"));
    BOOST_TEST_MESSAGE("FMM/LSM callable bond NPV " << npvFmm << ", stripped " << stripped << ", call value " << callFmm
                                                    << " +/- " << se);
    BOOST_CHECK(callFmm > 0.0);
    BOOST_CHECK_CLOSE(stripped - callFmm, npvFmm, 1e-8);

    // analytic straight bond from the trade's own cash flows: FMM's stripped NPV is the curve value
    Real straight = 0.0;
    for (const auto& cf : trade->legs().front())
        if (!cf->hasOccurred(kAsof))
            straight += cf->amount() * market->discountCurve("USD")->discount(cf->date());
    BOOST_TEST_MESSAGE("analytic straight bond " << straight);
    BOOST_CHECK_SMALL(stripped - straight, 1e-10);

    // ORE comparator: the LGM FD callable bond engine (24 time steps per year)
    auto lgmData = lgmEngineData("CallableBond");
    lgmData->engine("CallableBond") = "FD";
    lgmData->engineParameters("CallableBond") = {
        {"Scheme", "Douglas"}, {"StateGridPoints", "64"}, {"TimeStepsPerYear", "24"}, {"MesherEpsilon", "1e-4"}};
    auto lgmFactory = QuantLib::ext::make_shared<EngineFactory>(lgmData, market);
    auto lgmTrade = callableBondFromXml();
    lgmTrade->build(lgmFactory);
    const Real npvLgm = lgmTrade->instrument()->NPV();
    const Real callLgm = QuantLib::ext::any_cast<Real>(lgmTrade->instrument()->additionalResults().at("callPutValue"));
    const Real strippedLgm =
        QuantLib::ext::any_cast<Real>(lgmTrade->instrument()->additionalResults().at("strippedBondNpv"));
    BOOST_TEST_MESSAGE("LGM/FD callable bond NPV " << npvLgm << ", stripped " << strippedLgm << ", call value " << callLgm
                                                     << "; FMM - LGM call = " << callFmm - callLgm);
    BOOST_CHECK_SMALL(stripped - strippedLgm, 5e-5); // FD discretisation

    // ORE's event-time-grid variant (Grid builder, no intermediate time steps) omits the coupon paid at
    // the first event after the valuation date: at its own pay date the coupon has couponRatio 0 and is
    // excluded from the underlying, and no earlier time step revisits it (reported, not asserted)
    auto gridFactory = QuantLib::ext::make_shared<EngineFactory>(lgmEngineData("CallableBond"), market);
    auto gridTrade = callableBondFromXml();
    gridTrade->build(gridFactory);
    const Real strippedGrid =
        QuantLib::ext::any_cast<Real>(gridTrade->instrument()->additionalResults().at("strippedBondNpv"));
    const auto& firstCoupon = trade->legs().front().front();
    BOOST_TEST_MESSAGE("LGM/Grid (event-time grid) stripped bond NPV "
                       << strippedGrid << ": difference to FD " << strippedGrid - strippedLgm
                       << ", first coupon discounted "
                       << firstCoupon->amount() * market->discountCurve("USD")->discount(firstCoupon->date()));
    BOOST_CHECK_MESSAGE(std::fabs(callFmm - callLgm) < 0.10 * callLgm + 3.0 * se,
                        "FMM vs LGM callable bond option: " << callFmm << " vs " << callLgm);

    // issuer spread through the market's security spread
    market->spreadQuote_->setValue(0.005);
    for (auto& mb : fmmFactory->modelBuilders())
        mb.second->recalibrate();
    const Real npvFmmSpread = trade->instrument()->NPV();
    const Real strippedSpread = QuantLib::ext::any_cast<Real>(trade->instrument()->additionalResults().at("strippedBondNpv"));
    for (auto& mb : lgmFactory->modelBuilders())
        mb.second->recalibrate();
    const Real npvLgmSpread = lgmTrade->instrument()->NPV();
    const Real strippedLgmSpread =
        QuantLib::ext::any_cast<Real>(lgmTrade->instrument()->additionalResults().at("strippedBondNpv"));
    BOOST_TEST_MESSAGE("50 bp issuer spread: FMM " << npvFmmSpread << " (stripped " << strippedSpread << ") vs LGM "
                                                   << npvLgmSpread << " (stripped " << strippedLgmSpread << ")");
    BOOST_CHECK(strippedSpread < stripped);
    BOOST_CHECK_SMALL(strippedSpread - strippedLgmSpread, 5e-5);
    BOOST_CHECK_MESSAGE(std::fabs(npvFmmSpread - npvLgmSpread) < 0.10 * callLgm + 3.0 * se,
                        "FMM vs LGM callable bond with spread: " << npvFmmSpread << " vs " << npvLgmSpread);
}

BOOST_AUTO_TEST_CASE(testFmmNoticePeriods) {
    BOOST_TEST_MESSAGE("Testing notice periods through the engine factory: a Bermudan swaption with a one-month "
                       "notice (decision at the notice date, swap entry at the exercise date) and a callable bond "
                       "with a one-month notice in its CallData...");
    Settings::instance().evaluationDate() = kAsof;
    auto market = QuantLib::ext::make_shared<FmmTestMarket>(kAsof, 0.03, 0.0080);

    // swaption: ORE passes the notice dates as exercise dates and the settlement dates alongside
    auto plainFactory = QuantLib::ext::make_shared<EngineFactory>(fmmEngineData("BermudanSwaption", false), market);
    auto plain = bermudanPayerSwaption(0.031);
    plain->build(plainFactory);
    const Real npvPlain = plain->instrument()->NPV();
    const Real sePlain =
        QuantLib::ext::any_cast<Real>(plain->instrument()->additionalResults().at("fmmLsmLowerBoundStdError"));
    auto noticeFactory = QuantLib::ext::make_shared<EngineFactory>(fmmEngineData("BermudanSwaption", false), market);
    auto notice = bermudanPayerSwaption(0.031, "1M");
    notice->build(noticeFactory);
    const Real npvNotice = notice->instrument()->NPV();
    const auto& add = notice->instrument()->additionalResults();
    const Real seNotice = QuantLib::ext::any_cast<Real>(add.at("fmmLsmLowerBoundStdError"));
    const auto settle = QuantLib::ext::any_cast<std::vector<Date>>(add.at("fmmSettlementDates"));
    const auto exercise = QuantLib::ext::any_cast<std::vector<Date>>(add.at("fmmExerciseDates"));
    BOOST_REQUIRE_EQUAL(settle.size(), Size(4));
    BOOST_REQUIRE_EQUAL(exercise.size(), Size(4));
    for (Size i = 0; i < 4; ++i) {
        BOOST_CHECK_EQUAL(settle[i], Date(19, September, 2027 + static_cast<Year>(i)));
        BOOST_CHECK_EQUAL(exercise[i], Date(19, August, 2027 + static_cast<Year>(i)));
    }
    const Real periodsPlain = QuantLib::ext::any_cast<Real>(plain->instrument()->additionalResults().at("fmmGridPeriods"));
    const Real periodsNotice = QuantLib::ext::any_cast<Real>(add.at("fmmGridPeriods"));
    BOOST_TEST_MESSAGE("Bermudan payer: no notice " << npvPlain << " +/- " << sePlain << " (" << periodsPlain
                                                    << " grid periods), one-month notice " << npvNotice << " +/- "
                                                    << seNotice << " (" << periodsNotice << " grid periods)");
    // deciding a month early with less information cannot be worth more (beyond noise)
    BOOST_CHECK_MESSAGE(npvNotice < npvPlain + 3.0 * std::sqrt(sePlain * sePlain + seNotice * seNotice),
                        "notice-period value must not exceed the no-notice value: " << npvNotice << " vs " << npvPlain);
    BOOST_CHECK(npvNotice > 0.5 * npvPlain);
    // the notice dates are contractual grid dates (four extra periods) and the calibration helpers
    // expiring on them carry schedules anchored on the notice anniversaries, whose dates beyond the
    // mapping tolerance are added as well
    auto noticeBuilder = QuantLib::ext::dynamic_pointer_cast<FmmBuilder>(noticeFactory->modelBuilders().begin()->second);
    BOOST_REQUIRE(noticeBuilder);
    std::ostringstream gridDates;
    for (const Date& d : noticeBuilder->grid()->dates())
        gridDates << d << " ";
    BOOST_TEST_MESSAGE("notice grid: " << gridDates.str());
    BOOST_CHECK(periodsNotice >= periodsPlain + 4.0);

    // callable bond: NoticePeriod in the CallData, NoticePeriod model parameter for the grid
    const std::string noticeXml = "<NoticePeriod>1M</NoticePeriod><NoticeCalendar>NullCalendar</NoticeCalendar>"
                                  "<NoticeConvention>Unadjusted</NoticeConvention>";
    auto bondPlain = callableBondFromXml();
    auto bondPlainFactory = QuantLib::ext::make_shared<EngineFactory>(fmmEngineData("CallableBond", false), market);
    bondPlain->build(bondPlainFactory);
    const Real bondNpvPlain = bondPlain->instrument()->NPV();
    const Real bondSePlain =
        QuantLib::ext::any_cast<Real>(bondPlain->instrument()->additionalResults().at("fmmLsmLowerBoundStdError"));
    auto bondNotice = callableBondFromXml(noticeXml);
    BOOST_CHECK_EQUAL(bondNotice->data().callData().noticePeriod(), "1M");
    auto bondNoticeFactory = QuantLib::ext::make_shared<EngineFactory>(
        fmmEngineData("CallableBond", false, {{"NoticePeriod", "1M"}, {"NoticeCalendar", "NullCalendar"}}), market);
    bondNotice->build(bondNoticeFactory);
    const Real bondNpvNotice = bondNotice->instrument()->NPV();
    const auto& badd = bondNotice->instrument()->additionalResults();
    const Real bondSeNotice = QuantLib::ext::any_cast<Real>(badd.at("fmmLsmLowerBoundStdError"));
    const auto noticeDates = QuantLib::ext::any_cast<std::vector<Date>>(badd.at("fmmNoticeDates"));
    BOOST_REQUIRE_EQUAL(noticeDates.size(), Size(4));
    for (Size i = 0; i < 4; ++i)
        BOOST_CHECK_EQUAL(noticeDates[i], Date(19, August, 2027 + static_cast<Year>(i)));
    BOOST_TEST_MESSAGE("callable bond (holder): no notice " << bondNpvPlain << " +/- " << bondSePlain
                                                           << ", one-month notice " << bondNpvNotice << " +/- "
                                                           << bondSeNotice << ", call value "
                                                           << QuantLib::ext::any_cast<Real>(badd.at("callPutValue")));
    // the issuer's call with less information is worth less, so the holder's bond is worth more
    BOOST_CHECK_MESSAGE(bondNpvNotice >
                            bondNpvPlain - 3.0 * std::sqrt(bondSePlain * bondSePlain + bondSeNotice * bondSeNotice),
                        "notice-period bond value must not fall below the no-notice value: " << bondNpvNotice << " vs "
                                                                                             << bondNpvPlain);
    // the trade's notice dates must be grid dates: without the model parameter the engine refuses
    // (at pricing time, when the call data is mapped onto the grid)
    auto bondNoGrid = callableBondFromXml(noticeXml);
    auto bondNoGridFactory = QuantLib::ext::make_shared<EngineFactory>(fmmEngineData("CallableBond", false), market);
    bondNoGrid->build(bondNoGridFactory);
    BOOST_CHECK_THROW(bondNoGrid->instrument()->NPV(), std::exception);
    // XML round trip of the notice fields
    XMLDocument doc;
    XMLNode* node = bondNotice->toXML(doc);
    auto again = QuantLib::ext::make_shared<ore::data::CallableBond>();
    again->fromXML(node);
    BOOST_CHECK_EQUAL(again->data().callData().noticePeriod(), "1M");
    BOOST_CHECK_EQUAL(again->data().callData().noticeCalendar(), "NullCalendar");
    BOOST_CHECK_EQUAL(again->data().callData().noticeConvention(), "Unadjusted");
}

BOOST_AUTO_TEST_CASE(testFmmBuilderEvaluationDateMove) {
    BOOST_TEST_MESSAGE("Testing FmmBuilder on a moved evaluation date: the builder rebuilds the grid, the basket and "
                       "the calibration on the new reference date (engines built before the move must be rebuilt, "
                       "as for LGM)...");
    Settings::instance().evaluationDate() = kAsof;
    auto market = QuantLib::ext::make_shared<FmmTestMarket>(kAsof, 0.03, 0.0080);
    auto factory = QuantLib::ext::make_shared<EngineFactory>(fmmEngineData("BermudanSwaption", false), market);
    auto trade = bermudanPayerSwaption(0.031);
    trade->build(factory);
    trade->instrument()->NPV();
    auto builder = QuantLib::ext::dynamic_pointer_cast<FmmBuilder>(factory->modelBuilders().begin()->second);
    BOOST_REQUIRE(builder);
    const auto gridDates0 = builder->grid()->dates();
    BOOST_CHECK_EQUAL(builder->grid()->referenceDate(), kAsof);
    const Date moved = kAsof + 1 * Months;
    Settings::instance().evaluationDate() = moved; // flat curve and vol have settlement days 0: they follow
    BOOST_CHECK(builder->requiresRecalibration());
    builder->recalibrate();
    BOOST_CHECK_EQUAL(builder->grid()->referenceDate(), moved);
    BOOST_CHECK(builder->grid()->dates().front() == moved);
    BOOST_CHECK(builder->grid()->dates()[1] > moved);
    BOOST_CHECK(builder->grid()->dates().back() == gridDates0.back()); // same maturity anchor
    BOOST_CHECK(builder->grid()->numberOfRates() <= gridDates0.size() - 1); // lattice re-anchored, no new dates
    const auto& info = builder->calibrationInfo();
    BOOST_REQUIRE_EQUAL(info.modelVols.size(), Size(4));
    for (Size i = 0; i < info.modelVols.size(); ++i) {
        BOOST_CHECK_SMALL(info.modelVols[i] - 0.0080, 1e-7);
        BOOST_CHECK_CLOSE(info.expiryTimes[i],
                          Actual365Fixed().yearFraction(moved, Date(19, September, 2027 + static_cast<Year>(i))),
                          1e-6);
    }
    BOOST_TEST_MESSAGE("moved " << kAsof << " -> " << moved << ": grid " << gridDates0.size() - 1 << " -> "
                                << builder->grid()->numberOfRates() << " periods, first expiry time "
                                << info.expiryTimes.front());
    // the original trade is now seasoned (its swap accrues over the new valuation date), which the
    // FMM engines reject explicitly; a trade starting after the move prices on the new model
    auto seasonedFactory = QuantLib::ext::make_shared<EngineFactory>(fmmEngineData("BermudanSwaption", false), market);
    auto seasoned = bermudanPayerSwaption(0.031);
    seasoned->build(seasonedFactory); // the leg mapper runs at pricing time
    BOOST_CHECK_THROW(seasoned->instrument()->NPV(), std::exception);
    auto factory2 = QuantLib::ext::make_shared<EngineFactory>(fmmEngineData("BermudanSwaption", false), market);
    auto trade2 = bermudanPayerSwaption(0.031, "", "2026-12-19", "2031-12-19",
                                        {"2027-12-19", "2028-12-19", "2029-12-19", "2030-12-19"});
    trade2->build(factory2);
    const Real npv2 = trade2->instrument()->NPV();
    BOOST_TEST_MESSAGE("NPV of a forward-starting swaption one month later " << npv2);
    BOOST_CHECK(npv2 > 0.0);
}

BOOST_AUTO_TEST_CASE(testFmmJointCapFloorSwaptionCalibration) {
    BOOST_TEST_MESSAGE("Testing the cap/floor (optionlet) basket through the engine factory: joint bootstrap of the "
                       "ATM caplets from the market's optionlet surface and the coterminal swaptions...");
    Settings::instance().evaluationDate() = kAsof;
    auto market = QuantLib::ext::make_shared<FmmTestMarket>(kAsof, 0.03, 0.0080);

    // one factor, perfectly correlated, flat 80 bp caplets and 80 bp coterminal swaptions: the
    // swaptions are matched exactly (they are fitted last in each joint iteration) while the
    // backward-looking caplets keep a small structural residual, since their in-period variance is
    // not seen by the swaptions; the iteration reaches its cap and the residual is reported, not forced
    auto factory1 = QuantLib::ext::make_shared<EngineFactory>(
        fmmEngineData("BermudanSwaption", false,
                      {{"CapFloorBasket", "ATM"}, {"Factors", "1"}, {"RhoInf", "1.0"}, {"Beta", "0.0"},
                       {"JointMaxIterations", "15"}, {"JointToleranceBp", "0.05"}}),
        market);
    auto trade1 = bermudanPayerSwaption(0.031);
    trade1->build(factory1);
    const Real npv1 = trade1->instrument()->NPV();
    const auto& add1 = trade1->instrument()->additionalResults();
    BOOST_REQUIRE(add1.count("fmmCalibrationCapletMarketVols") > 0);
    const auto capMkt = QuantLib::ext::any_cast<std::vector<Real>>(add1.at("fmmCalibrationCapletMarketVols"));
    const auto capMdl = QuantLib::ext::any_cast<std::vector<Real>>(add1.at("fmmCalibrationCapletModelVols"));
    const auto capT = QuantLib::ext::any_cast<std::vector<Real>>(add1.at("fmmCalibrationCapletTimes"));
    const auto swpMkt = QuantLib::ext::any_cast<std::vector<Real>>(add1.at("fmmCalibrationMarketVols"));
    const auto swpMdl = QuantLib::ext::any_cast<std::vector<Real>>(add1.at("fmmCalibrationModelVols"));
    BOOST_REQUIRE_EQUAL(capMkt.size(), capMdl.size());
    BOOST_CHECK_EQUAL(capMkt.size(), Size(20)); // one caplet per quarterly period over the 5y grid
    BOOST_REQUIRE_EQUAL(swpMkt.size(), Size(4));
    Real worstCap = 0.0, worstSwp = 0.0;
    for (Size i = 0; i < capMkt.size(); ++i) {
        BOOST_CHECK_SMALL(capMkt[i] - 0.0080, 1e-12);
        worstCap = std::max(worstCap, std::fabs(capMdl[i] - capMkt[i]) * 1e4);
    }
    for (Size i = 0; i < swpMkt.size(); ++i)
        worstSwp = std::max(worstSwp, std::fabs(swpMdl[i] - swpMkt[i]) * 1e4);
    BOOST_TEST_MESSAGE("1F rhoInf=1: joint iterations " << QuantLib::ext::any_cast<Real>(add1.at("fmmJointIterations"))
                                                        << ", converged "
                                                        << QuantLib::ext::any_cast<bool>(add1.at("fmmJointConverged"))
                                                        << ", worst caplet residual " << worstCap
                                                        << " bp, worst swaption residual " << worstSwp << " bp, NPV "
                                                        << npv1 << ", first caplet time " << capT.front());
    // the alternation stops at its fixed point (stationary) before the 15-iteration cap; the
    // residual criterion is not met (structural, about 2 bp on this flat market)
    BOOST_CHECK(QuantLib::ext::any_cast<Real>(add1.at("fmmJointIterations")) <= 15.0);
    BOOST_CHECK(!QuantLib::ext::any_cast<bool>(add1.at("fmmJointConverged")));
    BOOST_CHECK(QuantLib::ext::any_cast<bool>(add1.at("fmmJointStationary")) ||
                QuantLib::ext::any_cast<Real>(add1.at("fmmJointIterations")) == 15.0);
    BOOST_TEST_MESSAGE("stationary " << QuantLib::ext::any_cast<bool>(add1.at("fmmJointStationary")) << ", last parameter change "
                                     << QuantLib::ext::any_cast<Real>(add1.at("fmmJointLastParameterChange")));
    BOOST_CHECK(worstSwp < 0.01);
    BOOST_CHECK(worstCap > 0.5 && worstCap < 3.0);
    BOOST_CHECK_CLOSE(QuantLib::ext::any_cast<Real>(add1.at("fmmCapletMaxResidualBp")), worstCap, 1e-6);
    BOOST_CHECK(npv1 > 0.0);

    // three factors with the default decorrelation: flat 80 / 80 targets are no longer consistent,
    // the iteration stops at its cap and the irreducible residual is reported, not forced
    auto factory3 = QuantLib::ext::make_shared<EngineFactory>(
        fmmEngineData("BermudanSwaption", false, {{"CapFloorBasket", "ATM"}, {"JointMaxIterations", "6"}, {"Tolerance", "0.05"}}),
        market);
    auto trade3 = bermudanPayerSwaption(0.031);
    trade3->build(factory3);
    const Real npv3 = trade3->instrument()->NPV();
    const auto& add3 = trade3->instrument()->additionalResults();
    const Real resid3 = QuantLib::ext::any_cast<Real>(add3.at("fmmCapletMaxResidualBp"));
    BOOST_TEST_MESSAGE("3F rhoInf=0.6: joint iterations " << QuantLib::ext::any_cast<Real>(add3.at("fmmJointIterations"))
                                                          << ", converged "
                                                          << QuantLib::ext::any_cast<bool>(add3.at("fmmJointConverged"))
                                                          << ", worst caplet residual " << resid3 << " bp, NPV " << npv3);
    BOOST_CHECK(resid3 >= 0.0);
    BOOST_CHECK(npv3 > 0.0);

    // XML round trip of the new data fields
    XMLDocument doc;
    doc.fromXMLString("<FMM key=\"USD-SOFR\"><CalibrationType>Bootstrap</CalibrationType><CapFloorBasket>ATM</CapFloorBasket>"
                      "<CapFloorHorizon>3Y</CapFloorHorizon><JointMaxIterations>7</JointMaxIterations>"
                      "<JointToleranceBp>0.2</JointToleranceBp><Volatility><InitialValue>0.003</InitialValue></Volatility>"
                      "<CalibrationSwaptions><Expiries>1Y</Expiries><Terms>4Y</Terms></CalibrationSwaptions></FMM>");
    FmmData data;
    data.fromXML(doc.getFirstNode("FMM"));
    BOOST_CHECK_EQUAL(data.capFloorBasket(), "ATM");
    BOOST_CHECK(data.capFloorHorizon() == 3 * Years);
    BOOST_CHECK_EQUAL(data.jointMaxIterations(), Size(7));
    BOOST_CHECK_CLOSE(data.jointToleranceBp(), 0.2, 1e-12);
    XMLDocument out;
    FmmData again;
    again.fromXML(data.toXML(out));
    BOOST_CHECK_EQUAL(again.capFloorBasket(), "ATM");
    BOOST_CHECK(again.capFloorHorizon() == 3 * Years);
    FmmData bad = data;
    bad.capFloorBasket() = "DealStrike";
    BOOST_CHECK_THROW(bad.validate(), QuantLib::Error);
}

BOOST_AUTO_TEST_CASE(testFmmPolicyModes) {
    BOOST_TEST_MESSAGE("Testing the exercise-policy treatment across revaluations (A6): retrained policy vs the "
                       "policy of the first valuation kept frozen, under a vol bump with recalibration...");
    Settings::instance().evaluationDate() = kAsof;
    auto market = QuantLib::ext::make_shared<FmmTestMarket>(kAsof, 0.03, 0.0080);
    auto retrainData = fmmEngineData("BermudanSwaption", false);
    auto frozenData = fmmEngineData("BermudanSwaption", false);
    frozenData->engineParameters("BermudanSwaption")["PolicyMode"] = "Frozen";
    auto retrainFactory = QuantLib::ext::make_shared<EngineFactory>(retrainData, market);
    auto frozenFactory = QuantLib::ext::make_shared<EngineFactory>(frozenData, market);
    auto retrain = bermudanPayerSwaption(0.031);
    auto frozen = bermudanPayerSwaption(0.031);
    retrain->build(retrainFactory);
    frozen->build(frozenFactory);
    const Real base = retrain->instrument()->NPV();
    const Real baseFrozen = frozen->instrument()->NPV();
    BOOST_CHECK_SMALL(base - baseFrozen, 1e-12); // identical first valuation (same seeds)
    BOOST_CHECK_EQUAL(QuantLib::ext::any_cast<std::string>(frozen->instrument()->additionalResults().at("fmmPolicyMode")),
                      "Frozen (policy trained)");
    BOOST_CHECK_EQUAL(QuantLib::ext::any_cast<std::string>(retrain->instrument()->additionalResults().at("fmmPolicyMode")),
                      "Retrain");

    market->volQuote_->setValue(0.0081); // 1 bp normal vol bump, recalibrated in both runs
    for (auto& mb : retrainFactory->modelBuilders())
        mb.second->recalibrate();
    for (auto& mb : frozenFactory->modelBuilders())
        mb.second->recalibrate();
    const Real bumped = retrain->instrument()->NPV();
    const Real bumpedFrozen = frozen->instrument()->NPV();
    const Real se = QuantLib::ext::any_cast<Real>(retrain->instrument()->additionalResults().at("fmmLsmLowerBoundStdError"));
    BOOST_CHECK_EQUAL(QuantLib::ext::any_cast<std::string>(frozen->instrument()->additionalResults().at("fmmPolicyMode")),
                      "Frozen (policy of the first valuation reused)");
    BOOST_TEST_MESSAGE("vega (1 bp, paired paths): retrained policy " << bumped - base << ", frozen policy "
                                                                      << bumpedFrozen - baseFrozen << " (LSM s.e. " << se
                                                                      << ")");
    // both vegas are positive and of the same order; the frozen policy is sub-optimal on the bumped
    // market, so its value cannot exceed the retrained one beyond the (paired, small) noise
    BOOST_CHECK(bumped > base);
    BOOST_CHECK(bumpedFrozen > baseFrozen);
    BOOST_CHECK(bumpedFrozen <= bumped + 0.1 * se);
    BOOST_CHECK(std::fabs((bumpedFrozen - baseFrozen) - (bumped - base)) < 0.5 * se);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
