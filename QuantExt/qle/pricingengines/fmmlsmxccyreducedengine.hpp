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

/*! \file qle/pricingengines/fmmlsmxccyreducedengine.hpp
    \brief Hagan-reduced pricing of a cancellable cross-currency swap on a single-currency FMM:
           the MultiLegOption (two currencies, rebates) is reduced to an Enter-style instrument in
           the modelled currency B (qle/models/fmmxccyreduction.hpp) and priced with the LSM
           policy and dual bound of the single-currency engine. The model is calibrated by the
           caller through a callback that receives the reduction's diagnostics (the strikes
           R_fix - s_k^eq and, for the improved method, the effective volatilities of Hagan
           eq. 3.41), so the trade-level calibration stays outside QuantExt's market-free scope.
*/

#pragma once

#include <qle/models/fmmxccyreduction.hpp>
#include <qle/pricingengines/fmmlsmmultilegoptionengine.hpp>

#include <functional>
#include <map>
#include <set>
#include <string>

namespace QuantExt {

//! what the caller's calibration returns: the calibrated B-currency model and the per-right
//! calibration quantities (in basket order = the reduction's rights)
struct FmmXccyCalibrationRecord {
    QuantLib::ext::shared_ptr<ForwardMarketModel> model;
    std::string method; //!< Standard | Improved
    std::vector<Real> strikes, marketVolsB, targetVols, sigmaA, sigmaFx, rhoFxB, rhoAB, omittedEffectVol;
    std::map<std::string, QuantLib::ext::any> extraResults;
};

class FmmLsmXccyReducedEngine : public MultiLegOption::engine {
public:
    //! the model grid on the trade's contractual dates
    using GridProvider = std::function<QuantLib::ext::shared_ptr<FmmGrid>(const std::set<Date>&)>;
    //! the calibrated model for the reduction's diagnostics on that grid
    using ModelProvider =
        std::function<FmmXccyCalibrationRecord(const FmmXccyReductionDiagnostics&, const FmmGrid&)>;
    FmmLsmXccyReducedEngine(const FmmXccyReductionInputs& inputs, const Currency& npvCurrency,
                            const GridProvider& gridProvider, const ModelProvider& modelProvider,
                            const FmmLsmEngineConfig& config);
    void calculate() const override;

private:
    FmmXccyReductionInputs inputs_;
    Currency npvCurrency_;
    GridProvider gridProvider_;
    ModelProvider modelProvider_;
    FmmLsmEngineConfig config_;
    mutable QuantLib::ext::optional<FmmLsmPolicy> frozenPolicy_;
};

} // namespace QuantExt
