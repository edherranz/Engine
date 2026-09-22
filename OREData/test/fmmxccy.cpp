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

#include <boost/test/unit_test.hpp>
#include <ored/utilities/toplevelfixture.hpp>

#include <ored/marketdata/marketimpl.hpp>
#include <ored/portfolio/builders/fmm.hpp>
#include <ored/portfolio/enginedata.hpp>
#include <ored/portfolio/enginefactory.hpp>
#include <ored/portfolio/swaption.hpp>
#include <ored/utilities/indexparser.hpp>
#include <ored/utilities/parsers.hpp>
#include <ored/utilities/xmlutils.hpp>

#include <qle/termstructures/flatcorrelation.hpp>

#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/termstructures/volatility/swaption/swaptionconstantvol.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include <cmath>
#include <map>
#include <sstream>

using namespace QuantLib;
using namespace QuantExt;
using namespace ore::data;
using std::make_pair;
using std::make_tuple;

namespace {

/* Two-currency test market of the cross-currency milestone: USD flat 4% (USD-SOFR), EUR flat
   2.5% (EUR-ESTER), EUR/USD spot 1.03 USD per EUR, flat FX vol 8%, constant normal swaption vols
   (USD 100 bp, EUR 75 bp), correlations in ORE's convention (FX-GENERIC-EUR-USD = USD per EUR):
   rate/rate 0.25, FX/EUR rate -0.25 (Hagan's rho(F^{A/B}, R^B) = +0.25), FX/USD rate 0. */
class XccyTestMarket : public MarketImpl {
public:
    XccyTestMarket(const Date& asof, const Real fxSpotUsdPerEur = 1.03, const Real rhoFxEurOre = -0.25)
        : MarketImpl(false) {
        asof_ = asof;
        const std::string cfg = Market::defaultConfiguration;
        auto flat = [&](const Real r) {
            return Handle<YieldTermStructure>(
                QuantLib::ext::make_shared<FlatForward>(0, NullCalendar(), r, Actual365Fixed()));
        };
        Handle<YieldTermStructure> ytsUsd = flat(0.04), ytsEur = flat(0.025);
        Handle<IborIndex> sofr(parseIborIndex("USD-SOFR", ytsUsd));
        Handle<IborIndex> ester(parseIborIndex("EUR-ESTER", ytsEur));
        iborIndices_[make_pair(cfg, "USD-SOFR")] = sofr;
        iborIndices_[make_pair(cfg, "EUR-ESTER")] = ester;
        for (Date d = asof - 10; d < asof; ++d) {
            if (sofr->isValidFixingDate(d))
                sofr->addFixing(d, 0.04);
            if (ester->isValidFixingDate(d))
                ester->addFixing(d, 0.025);
        }
        auto conventions = QuantLib::ext::make_shared<Conventions>();
        conventions->add(QuantLib::ext::make_shared<OisConvention>("USD-SOFR-OIS", "0", "USD-SOFR", "A365", "NullCalendar",
                                                                   "0", "false", "Annual", "F", "F", "Forward"));
        conventions->add(QuantLib::ext::make_shared<OisConvention>("EUR-ESTER-OIS", "0", "EUR-ESTER", "A365", "NullCalendar",
                                                                   "0", "false", "Annual", "F", "F", "Forward"));
        // the currency's OIS convention names its rate factor in the correlation curves
        conventions->add(QuantLib::ext::make_shared<OisConvention>("EUR-OIS", "0", "EUR-ESTER", "A365", "NullCalendar", "0",
                                                                   "false", "Annual", "F", "F", "Forward"));
        conventions->add(QuantLib::ext::make_shared<SwapIndexConvention>("USD-CMS-2Y", "USD-SOFR-OIS"));
        conventions->add(QuantLib::ext::make_shared<SwapIndexConvention>("USD-CMS-30Y", "USD-SOFR-OIS"));
        conventions->add(QuantLib::ext::make_shared<SwapIndexConvention>("EUR-CMS-2Y", "EUR-ESTER-OIS"));
        conventions->add(QuantLib::ext::make_shared<SwapIndexConvention>("EUR-CMS-30Y", "EUR-ESTER-OIS"));
        InstrumentConventions::instance().setConventions(conventions);

        yieldCurves_[make_tuple(cfg, YieldCurveType::Discount, "USD")] = ytsUsd;
        yieldCurves_[make_tuple(cfg, YieldCurveType::Discount, "EUR")] = ytsEur;
        yieldCurves_[make_tuple(cfg, YieldCurveType::Yield, "USD-SOFR")] = ytsUsd;
        yieldCurves_[make_tuple(cfg, YieldCurveType::Yield, "EUR-ESTER")] = ytsEur;
        addSwapIndex("USD-CMS-2Y", "USD-SOFR", cfg);
        addSwapIndex("USD-CMS-30Y", "USD-SOFR", cfg);
        addSwapIndex("EUR-CMS-2Y", "EUR-ESTER", cfg);
        addSwapIndex("EUR-CMS-30Y", "EUR-ESTER", cfg);
        for (const std::string& k : {"USD-SOFR", "USD"})
            swaptionIndexBases_[make_pair(cfg, k)] = make_pair("USD-CMS-2Y", "USD-CMS-30Y");
        for (const std::string& k : {"EUR-ESTER", "EUR"})
            swaptionIndexBases_[make_pair(cfg, k)] = make_pair("EUR-CMS-2Y", "EUR-CMS-30Y");
        auto constVol = [&](const Real v) {
            return Handle<QuantLib::SwaptionVolatilityStructure>(
                QuantLib::ext::make_shared<QuantLib::ConstantSwaptionVolatility>(0, NullCalendar(), Following, v,
                                                                                 Actual365Fixed(), Normal, 0.0));
        };
        swaptionCurves_[make_pair(cfg, "USD-SOFR")] = constVol(0.0100);
        swaptionCurves_[make_pair(cfg, "USD")] = constVol(0.0100);
        swaptionCurves_[make_pair(cfg, "EUR-ESTER")] = constVol(0.0075);
        swaptionCurves_[make_pair(cfg, "EUR")] = constVol(0.0075);
        // FX spot (USD per EUR), vol and the correlation curves
        std::map<std::string, Handle<Quote>> quotes;
        quotes["EURUSD"] = Handle<Quote>(QuantLib::ext::make_shared<SimpleQuote>(fxSpotUsdPerEur));
        fx_ = QuantLib::ext::make_shared<FXTriangulation>(quotes);
        fxVols_[make_pair(cfg, "EURUSD")] = Handle<BlackVolTermStructure>(
            QuantLib::ext::make_shared<BlackConstantVol>(0, NullCalendar(), 0.08, Actual365Fixed()));
        auto corr = [&](const Real rho) {
            return Handle<QuantExt::CorrelationTermStructure>(
                QuantLib::ext::make_shared<QuantExt::FlatCorrelation>(0, NullCalendar(), rho, Actual365Fixed()));
        };
        correlationCurves_[make_tuple(cfg, "USD-SOFR", "EUR-ESTER")] = corr(0.25);
        correlationCurves_[make_tuple(cfg, "FX-GENERIC-EUR-USD", "EUR-ESTER")] = corr(rhoFxEurOre);
        correlationCurves_[make_tuple(cfg, "FX-GENERIC-EUR-USD", "USD-SOFR")] = corr(0.0);
    }
};

const Date kAsof(19, September, 2026);

QuantLib::ext::shared_ptr<EngineData> xccyEngineData(const std::string& method, const std::string& product) {
    auto ed = QuantLib::ext::make_shared<EngineData>();
    ed->model(product) = "FMM";
    ed->engine(product) = "LSM";
    ed->modelParameters(product) = {{"Calibration", "Bootstrap"},
                                    {"CalibrationStrategy", "CoterminalDealStrike"},
                                    {"XccyMethod", method},
                                    {"Tolerance", "0.01"},
                                    {"Volatility", "0.0025"},
                                    {"VolatilityType", "DisplacedDiffusion"},
                                    {"Shift", "1/tau"},
                                    {"RhoInf", "0.6"},
                                    {"Beta", "0.08"},
                                    {"Factors", "3"},
                                    {"ApproximationMethod", "EffectiveShift"},
                                    {"Grid", "3M"},
                                    {"GridToleranceDays", "3"}};
    ed->engineParameters(product) = {{"TrainingPaths", "8192"}, {"ValuationPaths", "8192"}, {"TrainingSeed", "42"},
                                     {"ValuationSeed", "424242"}, {"BasisOrder", "2"}, {"DualBound", "false"}};
    ed->globalParameters()["GenerateAdditionalResults"] = "true";
    ed->model("Swap") = "DiscountedCashflows";
    ed->engine("Swap") = "DiscountingSwapEngine";
    ed->model("CrossCurrencySwap") = "DiscountedCashflows";
    ed->engine("CrossCurrencySwap") = "DiscountingCrossCurrencySwapEngine";
    return ed;
}

std::string legXml(const std::string& type, const bool payer, const std::string& ccy, const std::string& notional,
                   const std::string& body, const std::string& tenor, const bool finalExchange) {
    return "<LegData><LegType>" + type + "</LegType><Payer>" + (payer ? "true" : "false") + "</Payer><Currency>" + ccy +
           "</Currency><Notionals><Notional>" + notional + "</Notional>" +
           (finalExchange ? "<Exchanges><NotionalInitialExchange>false</NotionalInitialExchange>"
                            "<NotionalFinalExchange>true</NotionalFinalExchange>"
                            "<NotionalAmortizingExchange>false</NotionalAmortizingExchange></Exchanges>"
                          : "") +
           "</Notionals><DayCounter>A365</DayCounter><PaymentConvention>U</PaymentConvention>" + body +
           "<ScheduleData><Rules><StartDate>2026-09-19</StartDate><EndDate>2031-09-19</EndDate><Tenor>" + tenor +
           "</Tenor><Calendar>NullCalendar</Calendar><Convention>U</Convention><TermConvention>U</TermConvention>"
           "<Rule>Forward</Rule></Rules></ScheduleData></LegData>";
}

const std::string kExercise = "<ExerciseDates><ExerciseDate>2027-09-19</ExerciseDate><ExerciseDate>2028-09-19</ExerciseDate>"
                              "<ExerciseDate>2029-09-19</ExerciseDate><ExerciseDate>2030-09-19</ExerciseDate></ExerciseDates>";

// the cancellation right of the template: receive USD-SOFR + spread on 10m with the final exchange,
// pay EUR fixed 2.25% on 9.7m with the final exchange, fees +10m USD / -9.7m EUR on the exercise date
QuantLib::ext::shared_ptr<ore::data::Swaption> xccyCancelRight(const Real spread) {
    const std::string xml =
        "<Trade id=\"XCCY_RIGHT\"><TradeType>Swaption</TradeType><Envelope><CounterParty>CP</CounterParty>"
        "<NettingSetId>NS</NettingSetId></Envelope><SwaptionData><OptionData><LongShort>Long</LongShort>"
        "<OptionType>Call</OptionType><Style>Bermudan</Style><Settlement>Physical</Settlement>"
        "<PayOffAtExpiry>false</PayOffAtExpiry><ExerciseFees><ExerciseFee currency=\"USD\">10000000</ExerciseFee>"
        "<ExerciseFee currency=\"EUR\">-9700000</ExerciseFee></ExerciseFees>" +
        kExercise + "</OptionData>" +
        legXml("Floating", false, "USD", "10000000",
               "<FloatingLegData><Index>USD-SOFR</Index><IsInArrears>true</IsInArrears><Spreads><Spread>" +
                   std::to_string(spread) + "</Spread></Spreads></FloatingLegData>",
               "3M", true) +
        legXml("Fixed", true, "EUR", "9700000", "<FixedLegData><Rates><Rate>0.0225</Rate></Rates></FixedLegData>", "1Y",
               true) +
        "</SwaptionData></Trade>";
    auto trade = QuantLib::ext::make_shared<ore::data::Swaption>();
    XMLDocument doc;
    doc.fromXMLString(xml);
    trade->fromXML(doc.getFirstNode("Trade"));
    trade->id() = "XCCY_RIGHT"; // the portfolio loader sets the id; the engine cache is keyed by it
    return trade;
}

// the zero-spread limit: an EUR Bermudan payer swaption, fixed 2.25% against EUR-ESTER quarterly
QuantLib::ext::shared_ptr<ore::data::Swaption> eurBermudanPayer() {
    const std::string xml =
        "<Trade id=\"EUR_BERM\"><TradeType>Swaption</TradeType><Envelope><CounterParty>CP</CounterParty>"
        "<NettingSetId>NS</NettingSetId></Envelope><SwaptionData><OptionData><LongShort>Long</LongShort>"
        "<OptionType>Call</OptionType><Style>Bermudan</Style><Settlement>Physical</Settlement>"
        "<PayOffAtExpiry>false</PayOffAtExpiry>" +
        kExercise + "</OptionData>" +
        legXml("Fixed", true, "EUR", "9700000", "<FixedLegData><Rates><Rate>0.0225</Rate></Rates></FixedLegData>", "1Y",
               false) +
        legXml("Floating", false, "EUR", "9700000",
               "<FloatingLegData><Index>EUR-ESTER</Index><IsInArrears>true</IsInArrears><Spreads><Spread>0</Spread>"
               "</Spreads></FloatingLegData>",
               "3M", false) +
        "</SwaptionData></Trade>";
    auto trade = QuantLib::ext::make_shared<ore::data::Swaption>();
    XMLDocument doc;
    doc.fromXMLString(xml);
    trade->fromXML(doc.getFirstNode("Trade"));
    trade->id() = "EUR_BERM";
    return trade;
}

template <class T> T res(const QuantLib::ext::shared_ptr<ore::data::Trade>& t, const std::string& key) {
    const auto& add = t->instrument()->additionalResults();
    BOOST_REQUIRE_MESSAGE(add.count(key) > 0, "missing additional result " << key);
    return QuantLib::ext::any_cast<T>(add.at(key));
}

} // namespace

BOOST_FIXTURE_TEST_SUITE(OREDataTestSuite, ore::data::TopLevelFixture)

BOOST_AUTO_TEST_SUITE(FmmXccyTests)

BOOST_AUTO_TEST_CASE(testXccyCancelRightImprovedAndStandard) {
    BOOST_TEST_MESSAGE("A9.1: the cancellation right of a cross-currency swap through the engine factory (product "
                       "BermudanSwaption_XCcy, model FMM): Hagan's reduction on an EUR FMM, improved and standard "
                       "calibration targets, correlation sign convention, effective volatility formula...");
    Settings::instance().evaluationDate() = kAsof;
    auto market = QuantLib::ext::make_shared<XccyTestMarket>(kAsof);
    auto factory = QuantLib::ext::make_shared<EngineFactory>(xccyEngineData("Improved", "BermudanSwaption_XCcy"), market);
    auto trade = xccyCancelRight(0.0020);
    trade->build(factory);
    const Real npv = trade->instrument()->NPV();
    BOOST_TEST_MESSAGE("improved method NPV " << npv << " " << trade->npvCurrency());
    BOOST_CHECK_EQUAL(trade->npvCurrency(), "USD");
    BOOST_CHECK(npv > 0.0);
    BOOST_CHECK_EQUAL(res<std::string>(trade, "fmmXccyMethod"), "Improved");
    BOOST_CHECK_EQUAL(res<std::string>(trade, "fmmXccyCurrencyA"), "USD");
    BOOST_CHECK_EQUAL(res<std::string>(trade, "fmmXccyCurrencyB"), "EUR");
    BOOST_CHECK_EQUAL(res<std::string>(trade, "fmmXccyRateFactorB"), "EUR-ESTER");
    BOOST_CHECK_EQUAL(res<std::string>(trade, "fmmXccyFxIndexBperA"), "FX-GENERIC-USD-EUR");
    // ORE's fxRate() is the market spot moved to today with the two discount curves
    BOOST_CHECK_CLOSE(res<Real>(trade, "fmmXccyFxSpotBperA"), 1.0 / 1.03, 0.05);
    const auto strikes = res<std::vector<Real>>(trade, "fmmXccyStrikes");
    const auto sEq = res<std::vector<Real>>(trade, "fmmXccyEquivalentSpreads");
    const auto sigmaB = res<std::vector<Real>>(trade, "fmmXccyMarketVolsB");
    const auto target = res<std::vector<Real>>(trade, "fmmXccyTargetVols");
    const auto sigmaA = res<std::vector<Real>>(trade, "fmmXccyVolA");
    const auto sigmaFx = res<std::vector<Real>>(trade, "fmmXccyVolFx");
    const auto rhoFxB = res<std::vector<Real>>(trade, "fmmXccyRhoFxB");
    const auto rhoAB = res<std::vector<Real>>(trade, "fmmXccyRhoAB");
    const auto gammaA = res<std::vector<Real>>(trade, "fmmXccyGammaA");
    const auto gammaB = res<std::vector<Real>>(trade, "fmmXccyGammaB");
    const auto modelVols = res<std::vector<Real>>(trade, "fmmCalibrationModelVols");
    const auto marketVols = res<std::vector<Real>>(trade, "fmmCalibrationMarketVols");
    BOOST_REQUIRE_EQUAL(strikes.size(), Size(4));
    BOOST_REQUIRE_EQUAL(modelVols.size(), Size(4));
    for (Size k = 0; k < 4; ++k) {
        BOOST_TEST_MESSAGE("right " << k << ": s_eq " << sEq[k] * 1e4 << " bp, strike " << strikes[k] << ", sigmaB "
                                    << sigmaB[k] * 1e4 << " bp, target " << target[k] * 1e4 << " bp (sigmaA "
                                    << sigmaA[k] * 1e4 << " bp, sigmaFx " << sigmaFx[k] << ", rhoFxB " << rhoFxB[k]
                                    << ", rhoAB " << rhoAB[k] << ", gammaA " << gammaA[k] << ", gammaB " << gammaB[k]
                                    << "), calibrated model vol " << modelVols[k] * 1e4 << " bp");
        BOOST_CHECK(sEq[k] > 0.0 && sEq[k] < 0.0020 * 1.5);
        BOOST_CHECK_CLOSE(strikes[k], 0.0225 - sEq[k], 1e-8);
        BOOST_CHECK_CLOSE(sigmaB[k], 0.0075, 1e-10);
        BOOST_CHECK_CLOSE(sigmaA[k], 0.0100, 1e-10);
        BOOST_CHECK_CLOSE(sigmaFx[k], 0.08, 1e-10);
        // ORE quotes FX-GENERIC-EUR-USD (USD per EUR) against the EUR rate at -0.25; Hagan's
        // F^{A/B} is EUR per USD, so rho(F^{A/B}, R^B) = +0.25
        BOOST_CHECK_CLOSE(rhoFxB[k], 0.25, 1e-10);
        BOOST_CHECK_CLOSE(rhoAB[k], 0.25, 1e-10);
        // Hagan eq. 3.41
        const Real expected = sigmaB[k] + rhoFxB[k] * sigmaFx[k] * sEq[k] + gammaB[k] * sigmaB[k] * sEq[k] -
                              rhoAB[k] * gammaA[k] * sigmaA[k] * sEq[k];
        BOOST_CHECK_CLOSE(target[k], expected, 1e-8);
        BOOST_CHECK(target[k] > sigmaB[k]); // positive FX correlation dominates for these values
        // the calibration reproduces the effective vols
        BOOST_CHECK_SMALL(marketVols[k] - target[k], 1e-12);
        BOOST_CHECK_SMALL(modelVols[k] - target[k], 2e-6);
    }
    // standard method: the market vol at the strike is the target
    auto factoryStd = QuantLib::ext::make_shared<EngineFactory>(xccyEngineData("Standard", "BermudanSwaption_XCcy"), market);
    auto tradeStd = xccyCancelRight(0.0020);
    tradeStd->build(factoryStd);
    const Real npvStd = tradeStd->instrument()->NPV();
    const auto targetStd = res<std::vector<Real>>(tradeStd, "fmmXccyTargetVols");
    for (Size k = 0; k < 4; ++k)
        BOOST_CHECK_CLOSE(targetStd[k], 0.0075, 1e-10);
    const Real se = res<Real>(trade, "fmmLsmLowerBoundStdError") / (1.0 / 1.03);
    BOOST_TEST_MESSAGE("standard method NPV " << npvStd << " (improved - standard = " << npv - npvStd << ", s.e. " << se
                                              << ")");
    BOOST_CHECK(npv > npvStd); // higher target vols, higher option value
}

BOOST_AUTO_TEST_CASE(testXccyZeroSpreadEqualsEurBermudan) {
    BOOST_TEST_MESSAGE("A9.1: with zero spread the cancellation right priced by the reduced engine equals the EUR "
                       "Bermudan payer swaption priced by the single-currency FMM engine (same grid, calibration and "
                       "seeds), converted at the spot (Hagan section 3.1)...");
    Settings::instance().evaluationDate() = kAsof;
    auto market = QuantLib::ext::make_shared<XccyTestMarket>(kAsof);
    auto ed = xccyEngineData("Standard", "BermudanSwaption_XCcy");
    ed->model("BermudanSwaption") = "FMM";
    ed->engine("BermudanSwaption") = "LSM";
    ed->modelParameters("BermudanSwaption") = ed->modelParameters("BermudanSwaption_XCcy");
    ed->engineParameters("BermudanSwaption") = ed->engineParameters("BermudanSwaption_XCcy");
    auto factory = QuantLib::ext::make_shared<EngineFactory>(ed, market);
    auto xccy = xccyCancelRight(0.0);
    xccy->build(factory);
    auto eur = eurBermudanPayer();
    eur->build(factory);
    const Real npvXccy = xccy->instrument()->NPV(); // USD
    const Real npvEur = eur->instrument()->NPV();   // EUR
    const auto sEq = res<std::vector<Real>>(xccy, "fmmXccyEquivalentSpreads");
    for (const Real s : sEq)
        BOOST_CHECK_SMALL(s, 1e-14);
    BOOST_TEST_MESSAGE("reduced xccy right " << npvXccy << " USD = " << npvXccy / 1.03 << " EUR vs EUR Bermudan payer "
                                             << npvEur << " EUR (diff " << npvXccy / 1.03 - npvEur << ", s.e. "
                                             << res<Real>(eur, "fmmLsmLowerBoundStdError") << ")");
    auto show = [&](const std::string& what, const QuantLib::ext::shared_ptr<ore::data::Trade>& t) {
        const auto mv = res<std::vector<Real>>(t, "fmmCalibrationModelVols");
        const auto kv = res<std::vector<Real>>(t, "fmmCalibrationMarketVols");
        const auto pr = res<std::vector<Real>>(t, "fmmExerciseProbabilities");
        std::ostringstream o;
        for (Size k = 0; k < mv.size(); ++k)
            o << " " << kv[k] * 1e4 << "/" << mv[k] * 1e4;
        std::ostringstream q;
        for (const Real x : pr)
            q << " " << x;
        BOOST_TEST_MESSAGE(what << ": lower bound " << res<Real>(t, "fmmLsmLowerBound") << " +/- "
                                << res<Real>(t, "fmmLsmLowerBoundStdError") << ", underlying " << res<Real>(t, "underlyingNpv")
                                << ", grid periods " << res<Real>(t, "fmmGridPeriods") << ", max mismatch "
                                << res<Real>(t, "fmmGridMaxDateMismatchDays") << " days, market/model vols (bp)" << o.str()
                                << ", exercise probabilities" << q.str() << ", never "
                                << res<Real>(t, "fmmNoExerciseProbability"));
    };
    show("xccy reduced (EUR)", xccy);
    show("EUR Bermudan", eur);
    const auto fees = res<std::vector<Real>>(xccy, "fmmXccyFeeFlowsB");
    const auto strikes = res<std::vector<Real>>(xccy, "fmmXccyStrikes");
    for (Size k = 0; k < fees.size(); ++k)
        BOOST_TEST_MESSAGE("right " << k << ": fee " << fees[k] << " EUR, strike " << strikes[k]);
    BOOST_CHECK_EQUAL(res<Real>(xccy, "fmmGridPeriods"), res<Real>(eur, "fmmGridPeriods"));
    BOOST_CHECK_SMALL(res<Real>(xccy, "fmmXccyValueB") - npvEur, 1e-6 * 9700000.0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
