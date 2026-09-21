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

#include <orea/app/analytics/fmmexposureanalytic.hpp>

#include <orea/app/inputparameters.hpp>

#include <orea/cube/cube_io.hpp>
#include <orea/cube/inmemorycube.hpp>
#include <orea/scenario/aggregationscenariodata.hpp>

#include <ored/model/fmmbuilder.hpp>
#include <ored/portfolio/callablebond.hpp>
#include <ored/portfolio/portfolio.hpp>
#include <ored/portfolio/trade.hpp>
#include <ored/report/inmemoryreport.hpp>
#include <ored/utilities/log.hpp>
#include <ored/utilities/parsers.hpp>
#include <ored/utilities/to_string.hpp>

#include <qle/instruments/callablebond.hpp>
#include <qle/instruments/multilegoption.hpp>
#include <qle/models/fmmexposure.hpp>
#include <qle/pricingengines/fmmlsmcallablebondengine.hpp>
#include <qle/pricingengines/fmmlsmmultilegoptionengine.hpp>

#include <ql/cashflows/coupon.hpp>
#include <ql/instruments/swap.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

using namespace ore::data;
using namespace QuantLib;
using namespace QuantExt;

namespace ore {
namespace analytics {

void FmmExposureVariables::loadVariablesImpl(const QuantLib::ext::shared_ptr<InputParameters>& inputs) {
    inputs->loadParameterXML<ScenarioGeneratorData>(scenarioGeneratorData_, "simulation", "scenarioGeneratorData");
    if (!scenarioGeneratorData_)
        inputs->loadParameterXML<ScenarioGeneratorData>(scenarioGeneratorData_, "simulation", "simulationConfigFile",
                                                        true);
    inputs->loadParameterXML<ore::data::FmmData>(scenarioModel_, "fmmExposure", "fmmScenarioModelFile", true);
    auto sizeParser = [](const std::string& s) { return static_cast<Size>(parseInteger(s)); };
    auto seedParser = [](const std::string& s) { return static_cast<BigNatural>(parseInteger(s)); };
    inputs->loadParameter<Size>(trainingPaths_, "fmmExposure", "trainingPaths", false, sizeParser);
    inputs->loadParameter<BigNatural>(trainingSeed_, "fmmExposure", "trainingSeed", false, seedParser);
    inputs->loadParameter<std::string>(trainingSequence_, "fmmExposure", "trainingSequence", false);
    inputs->loadParameter<Size>(policyTrainingPaths_, "fmmExposure", "policyTrainingPaths", false, sizeParser);
    inputs->loadParameter<Size>(policyValuationPaths_, "fmmExposure", "policyValuationPaths", false, sizeParser);
    inputs->loadParameter<BigNatural>(policyTrainingSeed_, "fmmExposure", "policyTrainingSeed", false, seedParser);
    inputs->loadParameter<BigNatural>(policyValuationSeed_, "fmmExposure", "policyValuationSeed", false, seedParser);
    inputs->loadParameter<Size>(regressionOrder_, "fmmExposure", "regressionOrder", false, sizeParser);
    inputs->loadParameter<Size>(minRegressionPaths_, "fmmExposure", "minRegressionPaths", false, sizeParser);
    inputs->loadParameter<bool>(includeFlowsOnExposureDate_, "fmmExposure", "includeFlowsOnExposureDate", false,
                                [](const std::string& s) { return parseBool(s); });
    inputs->loadParameter<double>(calibrationTolerance_, "fmmExposure", "calibrationTolerance", false,
                                  [](const std::string& s) { return parseReal(s); });
}

void FmmExposureAnalyticImpl::setUpConfigurations() {
    auto vars = QuantLib::ext::dynamic_pointer_cast<FmmExposureVariables>(inputVariables_);
    analytic()->configurations().todaysMarketParams = inputs_->todaysMarketParams();
    analytic()->configurations().scenarioGeneratorData = vars->scenarioGeneratorData_;
}

namespace {

struct TradeInfo {
    std::string id;
    QuantLib::ext::shared_ptr<Trade> trade;
    enum class Kind { Swap, Swaption, CallableBond } kind = Kind::Swap;
    Real issuerSpread = 0.0;
    Real sign = 1.0; //!< applied to the engine's values for the cube (issuer-view structures)
    std::string perspective;
};

void addLegDates(const Leg& leg, const Date& asof, std::set<Date>& dates) {
    for (const auto& cf : leg) {
        if (cf->date() > asof)
            dates.insert(cf->date());
        if (auto c = QuantLib::ext::dynamic_pointer_cast<Coupon>(cf)) {
            if (c->accrualStartDate() > asof)
                dates.insert(c->accrualStartDate());
            if (c->accrualEndDate() > asof)
                dates.insert(c->accrualEndDate());
        }
    }
}

std::string anyToString(const QuantLib::ext::any& a) {
    if (a.type() == typeid(Real))
        return ore::data::to_string(QuantLib::ext::any_cast<Real>(a));
    if (a.type() == typeid(bool))
        return QuantLib::ext::any_cast<bool>(a) ? "true" : "false";
    if (a.type() == typeid(std::string))
        return QuantLib::ext::any_cast<std::string>(a);
    if (a.type() == typeid(Size))
        return std::to_string(QuantLib::ext::any_cast<Size>(a));
    if (a.type() == typeid(int))
        return std::to_string(QuantLib::ext::any_cast<int>(a));
    if (a.type() == typeid(Date))
        return ore::data::to_string(QuantLib::ext::any_cast<Date>(a));
    return "";
}

} // namespace

void FmmExposureAnalyticImpl::runAnalytic(const QuantLib::ext::shared_ptr<InMemoryLoader>& loader,
                                          const std::set<std::string>& runTypes) {
    auto vars = QuantLib::ext::dynamic_pointer_cast<FmmExposureVariables>(inputVariables_);
    QL_REQUIRE(vars && vars->scenarioGeneratorData_, "FmmExposureAnalytic: simulation configuration missing");
    QL_REQUIRE(vars->scenarioModel_, "FmmExposureAnalytic: fmmScenarioModelFile missing");
    LOG("FmmExposureAnalytic::runAnalytic called");
    CONSOLEW("FMM_EXPOSURE: Build Market");
    analytic()->buildMarket(loader);
    CONSOLE("OK");
    CONSOLEW("FMM_EXPOSURE: Build Portfolio");
    analytic()->buildPortfolio();
    CONSOLE("OK");

    const Date asof = inputs_->asof();
    const auto market = analytic()->market();
    const std::string marketConfig = inputs_->marketConfig("pricing");
    const auto& sgd = vars->scenarioGeneratorData_;
    QL_REQUIRE(!sgd->withCloseOutLag(),
               "FmmExposureAnalytic: a close-out lag is not supported in this release (baseline: depth 1 cube)");
    const auto grid = sgd->getGrid();
    const std::vector<Date>& dates = grid->valuationDates();
    QL_REQUIRE(!dates.empty(), "FmmExposureAnalytic: empty exposure date grid");
    const Size samples = sgd->samples();
    auto data = QuantLib::ext::make_shared<FmmData>(*vars->scenarioModel_);
    const Natural tol = data->gridToleranceDays();
    Actual365Fixed dc;

    // 1. the trades' instruments and contractual dates
    std::vector<TradeInfo> infos;
    std::set<Date> contractual;
    for (const auto& [id, trade] : analytic()->portfolio()->trades()) {
        auto ql = trade->instrument()->qlInstrument();
        TradeInfo ti;
        ti.id = id;
        ti.trade = trade;
        if (QuantLib::ext::dynamic_pointer_cast<QuantExt::MultiLegOption>(ql)) {
            ti.kind = TradeInfo::Kind::Swaption;
            ti.perspective = "option holder (Enter)";
            MultiLegOption::arguments args;
            ql->setupArguments(&args); // through the base class: some overrides are private
            for (const auto& leg : args.legs)
                addLegDates(leg, asof, contractual);
            if (args.exercise)
                for (const Date& d : args.exercise->dates())
                    if (d > asof)
                        contractual.insert(d);
            for (const Date& d : args.settlementDates)
                if (d > asof)
                    contractual.insert(d);
        } else if (QuantLib::ext::dynamic_pointer_cast<QuantExt::CallableBond>(ql)) {
            ti.kind = TradeInfo::Kind::CallableBond;
            ti.sign = -1.0; // the engine values the issuer's cancellable structure; the cube holds the holder's value
            ti.perspective = "bond holder (issuer Cancel structure negated)";
            QuantExt::CallableBond::arguments args;
            ql->setupArguments(&args);
            addLegDates(args.cashflows, asof, contractual);
            for (const auto& cd : args.callData) {
                if (cd.exerciseDate > asof)
                    contractual.insert(cd.exerciseDate);
                if (cd.noticeDate != Date() && cd.noticeDate > asof)
                    contractual.insert(cd.noticeDate);
            }
            if (auto oreBond = QuantLib::ext::dynamic_pointer_cast<ore::data::CallableBond>(trade)) {
                try {
                    ti.issuerSpread = market->securitySpread(oreBond->bondData().securityId(), marketConfig)->value();
                } catch (...) {
                    ti.issuerSpread = 0.0;
                }
            }
        } else if (QuantLib::ext::dynamic_pointer_cast<QuantLib::Swap>(ql)) {
            ti.kind = TradeInfo::Kind::Swap;
            ti.perspective = "swap holder";
            QuantLib::Swap::arguments args;
            ql->setupArguments(&args);
            for (const auto& leg : args.legs)
                addLegDates(leg, asof, contractual);
        } else {
            QL_FAIL("FmmExposureAnalytic: trade " << id << " (" << trade->tradeType()
                                                  << ") is not supported (swaps, physically settled Bermudan "
                                                     "swaptions and callable bonds)");
        }
        infos.push_back(ti);
    }
    QL_REQUIRE(!infos.empty(), "FmmExposureAnalytic: empty portfolio");
    contractual.insert(dates.back()); // the model grid covers the exposure horizon
    // contractual dates within the mapping tolerance of an earlier one collapse (no tiny periods;
    // the leg mapper absorbs the offset)
    std::vector<Date> merged;
    for (const Date& d : contractual)
        if (merged.empty() || d - merged.back() > static_cast<Integer>(tol))
            merged.push_back(d);
    data->gridDates() = merged;

    // 2. the scenario model of the netting set
    CONSOLEW("FMM_EXPOSURE: Calibrate scenario model");
    auto builder = QuantLib::ext::make_shared<FmmBuilder>(market, data, marketConfig, vars->calibrationTolerance_, false,
                                                          "", true, "FMM_EXPOSURE_SCENARIO_MODEL");
    auto model = builder->modelAsFmm();
    auto fmmGrid = builder->grid();
    CONSOLE("OK");
    LOG("FmmExposureAnalytic: scenario model on " << fmmGrid->numberOfRates() << " periods to "
                                                  << fmmGrid->dates().back() << ", calibration error "
                                                  << builder->error());

    // 3. instruments on the common grid and the exercise policies
    CONSOLEW("FMM_EXPOSURE: Map trades and train policies");
    FmmLsmConfig lsm;
    lsm.trainingPaths = vars->policyTrainingPaths_;
    lsm.valuationPaths = vars->policyValuationPaths_;
    lsm.trainingSeed = vars->policyTrainingSeed_;
    lsm.valuationSeed = vars->policyValuationSeed_;
    lsm.basisOrder = vars->regressionOrder_;
    std::vector<FmmExposureTrade> trades;
    std::vector<Real> lsmValue(infos.size(), Null<Real>()), lsmValueSe(infos.size(), Null<Real>());
    for (Size i = 0; i < infos.size(); ++i) {
        auto& ti = infos[i];
        auto ql = ti.trade->instrument()->qlInstrument();
        FmmExposureTrade t;
        t.id = ti.id;
        if (ti.kind == TradeInfo::Kind::Swaption) {
            MultiLegOption::arguments args;
            ql->setupArguments(&args);
            std::vector<Date> usedSettle;
            t.instrument = fmmMapMultiLegOption(args, *fmmGrid, tol, asof, usedSettle);
        } else if (ti.kind == TradeInfo::Kind::CallableBond) {
            QuantExt::CallableBond::arguments args;
            ql->setupArguments(&args);
            std::vector<Date> noticeDates;
            t.instrument = fmmMapCallableBond(args, *fmmGrid, tol, asof, ti.issuerSpread, noticeDates);
        } else {
            QuantLib::Swap::arguments args;
            ql->setupArguments(&args);
            t.instrument = fmmMapLegs(args.legs, args.payer, *fmmGrid, tol, asof);
        }
        if (!t.instrument.rights.empty()) {
            FmmLsmPricer pricer(model, t.instrument, lsm);
            const auto res = pricer.calculate();
            t.policy = pricer.policy();
            const bool cancel = t.instrument.style == FmmCallableInstrument::Style::Cancel;
            lsmValue[i] = ti.sign * (cancel ? res.lowerBoundCv : res.lowerBound);
            lsmValueSe[i] = cancel ? res.lowerBoundCvSe : res.lowerBoundSe;
        }
        trades.push_back(t);
    }
    CONSOLE("OK");

    // 4. exposure simulation
    CONSOLEW("FMM_EXPOSURE: Simulate exposure cube");
    std::vector<Time> times;
    for (const Date& d : dates)
        times.push_back(dc.yearFraction(asof, d));
    FmmExposureConfig cfg;
    cfg.trainingPaths = vars->trainingPaths_;
    cfg.trainingSeed = vars->trainingSeed_;
    cfg.trainingSequence = parseSequenceType(vars->trainingSequence_);
    cfg.exposurePaths = samples;
    cfg.exposureSeed = sgd->seed();
    cfg.exposureSequence = sgd->sequenceType();
    cfg.basisOrder = vars->regressionOrder_;
    cfg.minRegressionPaths = vars->minRegressionPaths_;
    cfg.includeFlowsOnExposureDate = vars->includeFlowsOnExposureDate_;
    FmmExposureEngine engine(model, trades, times, cfg);
    const auto res = engine.calculate();
    CONSOLE("OK");
    LOG("FmmExposureAnalytic: exposure simulation of " << trades.size() << " trades, " << dates.size() << " dates, "
                                                        << samples << " samples in " << res.runtimeSeconds << " s");

    // 5. cube and aggregation scenario data (native formats, written by the application)
    std::set<std::string> ids;
    for (const auto& ti : infos)
        ids.insert(ti.id);
    auto cube = QuantLib::ext::make_shared<DoublePrecisionInMemoryCube>(asof, ids, dates, samples, 1, 0.0);
    auto asd = QuantLib::ext::make_shared<InMemoryAggregationScenarioData>(dates.size(), samples);
    for (Size i = 0; i < infos.size(); ++i) {
        const Size idx = cube->idsAndIndexes().at(infos[i].id);
        cube->setT0(infos[i].sign * res.t0Values[i], idx, 0);
        for (Size d = 0; d < dates.size(); ++d)
            for (Size n = 0; n < samples; ++n)
                cube->set(infos[i].sign * res.deflatedValues[i][d * samples + n], idx, d, n, 0);
    }
    for (Size d = 0; d < dates.size(); ++d)
        for (Size n = 0; n < samples; ++n)
            asd->set(d, n, res.numeraire[d * samples + n], AggregationScenarioDataType::Numeraire, "");
    analytic()->npvCubes()[LABEL]["cube"] =
        QuantLib::ext::make_shared<NPVCubeWithMetaData>(cube, sgd, false, QuantLib::ext::nullopt);
    analytic()->mktCubes()[LABEL]["scenariodata"] = asd;

    // 6. reports
    auto t0Report = QuantLib::ext::make_shared<InMemoryReport>(inputs_->reportBufferSize());
    t0Report->addColumn("TradeId", std::string())
        .addColumn("TradeType", std::string())
        .addColumn("Perspective", std::string())
        .addColumn("Rights", Size())
        .addColumn("LastFlowDate", Date())
        .addColumn("PricingEngineNPV", double(), 4)
        .addColumn("ScenarioModelLsmValue", double(), 4)
        .addColumn("ScenarioModelLsmStdError", double(), 4)
        .addColumn("CubeT0", double(), 4)
        .addColumn("CubeT0StdError", double(), 4)
        .addColumn("CurveValue", double(), 4)
        .addColumn("IssuerSpread", double(), 8);
    for (Size i = 0; i < infos.size(); ++i) {
        Real pricingNpv = Null<Real>();
        try {
            pricingNpv = infos[i].trade->instrument()->NPV();
        } catch (const std::exception& e) {
            ALOG("FmmExposureAnalytic: pricing engine NPV of " << infos[i].id << " failed: " << e.what());
        }
        t0Report->next()
            .add(infos[i].id)
            .add(infos[i].trade->tradeType())
            .add(infos[i].perspective)
            .add(trades[i].instrument.rights.size())
            .add(fmmGrid->dates()[trades[i].instrument.lastFlowIdx])
            .add(pricingNpv)
            .add(lsmValue[i])
            .add(lsmValueSe[i])
            .add(infos[i].sign * res.t0Values[i])
            .add(res.t0ValuesSe[i])
            .add(infos[i].sign * res.t0CurveValues[i])
            .add(infos[i].issuerSpread);
    }
    t0Report->end();
    analytic()->addReport(LABEL, "fmmexposure_t0", t0Report);

    auto diagReport = QuantLib::ext::make_shared<InMemoryReport>(inputs_->reportBufferSize());
    diagReport->addColumn("TradeId", std::string())
        .addColumn("Date", Date())
        .addColumn("Time", double(), 6)
        .addColumn("AliveFraction", double(), 6)
        .addColumn("NoticedFraction", double(), 6)
        .addColumn("SettledFraction", double(), 6)
        .addColumn("ExercisedFraction", double(), 6)
        .addColumn("RegressionPaths", Size())
        .addColumn("RegressionR2", double(), 6)
        .addColumn("RegressionResidualSd", double(), 4)
        .addColumn("OutOfSampleBias", double(), 4)
        .addColumn("OutOfSampleBiasSe", double(), 4)
        .addColumn("MartingaleCheck", double(), 4)
        .addColumn("MartingaleCheckSe", double(), 4)
        .addColumn("T0Value", double(), 4);
    for (Size i = 0; i < infos.size(); ++i) {
        for (Size d = 0; d < dates.size(); ++d) {
            const auto& dg = res.diagnostics[i][d];
            diagReport->next()
                .add(infos[i].id)
                .add(dates[d])
                .add(times[d])
                .add(dg.aliveFraction)
                .add(dg.noticedFraction)
                .add(dg.settledFraction)
                .add(dg.exercisedFraction)
                .add(dg.regressionPaths)
                .add(dg.regressionR2)
                .add(dg.regressionResidualSd)
                .add(dg.outOfSampleBias)
                .add(dg.outOfSampleBiasSe)
                .add(dg.martingaleCheck)
                .add(dg.martingaleCheckSe)
                .add(res.t0Values[i]);
        }
    }
    diagReport->end();
    analytic()->addReport(LABEL, "fmmexposure_diagnostics", diagReport);

    auto calReport = QuantLib::ext::make_shared<InMemoryReport>(inputs_->reportBufferSize());
    calReport->addColumn("Key", std::string()).addColumn("Value", std::string());
    for (const auto& kv : builder->calibrationInfo().additionalResults()) {
        const std::string s = anyToString(kv.second);
        if (!s.empty()) {
            calReport->next().add(kv.first).add(s);
            continue;
        }
        if (kv.second.type() == typeid(std::vector<Real>)) {
            const auto& v = QuantLib::ext::any_cast<std::vector<Real>>(kv.second);
            for (Size k = 0; k < v.size(); ++k)
                calReport->next().add(kv.first + "[" + std::to_string(k) + "]").add(ore::data::to_string(v[k]));
        } else if (kv.second.type() == typeid(std::vector<std::string>)) {
            const auto& v = QuantLib::ext::any_cast<std::vector<std::string>>(kv.second);
            for (Size k = 0; k < v.size(); ++k)
                calReport->next().add(kv.first + "[" + std::to_string(k) + "]").add(v[k]);
        }
    }
    calReport->next().add("gridPeriods").add(std::to_string(fmmGrid->numberOfRates()));
    calReport->next().add("gridLastDate").add(ore::data::to_string(fmmGrid->dates().back()));
    calReport->next().add("exposureDates").add(std::to_string(dates.size()));
    calReport->next().add("samples").add(std::to_string(samples));
    calReport->next().add("runtimeSeconds").add(ore::data::to_string(res.runtimeSeconds));
    calReport->end();
    analytic()->addReport(LABEL, "fmmexposure_calibration", calReport);
}

} // namespace analytics
} // namespace ore
