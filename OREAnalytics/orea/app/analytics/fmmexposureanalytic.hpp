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

/*! \file orea/app/analytics/fmmexposureanalytic.hpp
    \brief FMM exposure analytic (A7 route 1): simulates a netting set on one calibrated FMM and
           writes the native NPV cube and aggregation scenario data for ORE's XVA post-processor
*/

#pragma once

#include <orea/app/analytic.hpp>
#include <orea/scenario/scenariogeneratordata.hpp>
#include <ored/model/fmmdata.hpp>

namespace ore {
namespace analytics {

struct FmmExposureVariables : public InputVariables {
    void loadVariablesImpl(const QuantLib::ext::shared_ptr<InputParameters>& inputs) override;

    QuantLib::ext::shared_ptr<ScenarioGeneratorData> scenarioGeneratorData_;
    QuantLib::ext::shared_ptr<ore::data::FmmData> scenarioModel_;
    std::string scenarioModelFile_;
    QuantLib::Size trainingPaths_ = 16384;
    QuantLib::BigNatural trainingSeed_ = 42;
    std::string trainingSequence_ = "MersenneTwister";
    QuantLib::Size policyTrainingPaths_ = 16384;
    QuantLib::Size policyValuationPaths_ = 16384;
    QuantLib::BigNatural policyTrainingSeed_ = 42;
    QuantLib::BigNatural policyValuationSeed_ = 4242;
    QuantLib::Size regressionOrder_ = 2;
    QuantLib::Size minRegressionPaths_ = 100;
    bool includeFlowsOnExposureDate_ = false;
    double calibrationTolerance_ = 0.001;
};

/*! ore.xml section "fmmExposure" (active Y): parameters fmmScenarioModelFile (FmmData XML, root
    <FMM key="...">: the scenario model of the netting set with its calibration basket), trainingPaths,
    trainingSeed, trainingSequence (conditional-value regressions), policyTrainingPaths,
    policyValuationPaths, policyTrainingSeed, policyValuationSeed (A4 exercise policies),
    regressionOrder, minRegressionPaths, includeFlowsOnExposureDate, calibrationTolerance. The
    exposure grid, samples, sequence and seed come from the simulation configuration
    (simulation/simulationConfigFile), the file names from simulation/cubeFile and
    simulation/aggregationScenarioDataFileName. The analytic prices the portfolio with the
    configured pricing engines (per-trade pricing calibrations, reported next to the scenario
    model's values), builds one FMM on the union grid of the netting set, trains the policies and
    the conditional-value estimators, simulates the exposure paths and registers the cube (depth 1,
    values V(t) / B(t) in base currency) and the aggregation data (Numeraire = B(t)). Reports:
    fmmexposure_t0, fmmexposure_diagnostics, fmmexposure_calibration. */
class FmmExposureAnalyticImpl : public Analytic::Impl {
public:
    static constexpr const char* LABEL = "FMM_EXPOSURE";
    explicit FmmExposureAnalyticImpl(const QuantLib::ext::shared_ptr<InputParameters>& inputs)
        : Analytic::Impl(inputs, QuantLib::ext::make_shared<FmmExposureVariables>()) {
        setLabel(LABEL);
    }
    void runAnalytic(const QuantLib::ext::shared_ptr<ore::data::InMemoryLoader>& loader,
                     const std::set<std::string>& runTypes = {}) override;
    void setUpConfigurations() override;
};

class FmmExposureAnalytic : public Analytic {
public:
    FmmExposureAnalytic(const QuantLib::ext::shared_ptr<InputParameters>& inputs,
                        const QuantLib::ext::weak_ptr<ore::analytics::AnalyticsManager>& analyticsManager)
        : Analytic(std::make_unique<FmmExposureAnalyticImpl>(inputs), {"FMM_EXPOSURE"}, inputs, analyticsManager, false,
                   false, false, false) {}
};

} // namespace analytics
} // namespace ore
