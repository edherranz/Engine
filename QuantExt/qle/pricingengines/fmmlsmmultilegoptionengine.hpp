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

/*! \file fmmlsmmultilegoptionengine.hpp
    \brief MultiLegOption engine on the FMM Longstaff-Schwartz pricer: Bermudan swaptions and
           cancellable swaps (as swap + option to enter the offsetting legs) on compounded-RFR or
           grid-period ibor legs, physically settled, grid-aligned exercise (ore-fmm
           docs/A5_PLAN.md sections 3-5). Unsupported conventions are rejected explicitly.
    \ingroup engines
*/

#ifndef quantext_fmm_lsm_multileg_option_engine_hpp
#define quantext_fmm_lsm_multileg_option_engine_hpp

#include <qle/instruments/multilegoption.hpp>
#include <qle/models/fmmgrid.hpp>
#include <qle/models/fmmlsmpricer.hpp>

namespace QuantExt {

//! LSM engine configuration shared by the FMM trade engines
struct FmmLsmEngineConfig {
    FmmLsmConfig lsm;
    bool dualBound = false; //!< also run the Andersen-Broadie dual bound (nested simulation)
    Size dualOuterPaths = 512, dualInnerPaths = 64;
    BigNatural dualSeed = 20260920;
    Natural gridToleranceDays = 3;
};

//! maps ORE / QuantLib legs onto an FmmCallableInstrument flow set (shared by the trade engines)
class FmmLegMapper {
public:
    FmmLegMapper(const FmmGrid& grid, const Natural toleranceDays, const Date& today);
    //! add a leg with the given sign (+1 receive, -1 pay) from the instrument owner's view
    void addLeg(const Leg& leg, const Real sign, FmmCallableInstrument& inst, const std::string& what) const;
    //! t = 0 curve value of the instrument's flows (fixed, grid floats, compounded floats)
    static Real underlyingCurveValue(const FmmCallableInstrument& inst, const FmmParametrization& p);

private:
    const FmmGrid& grid_;
    Natural toleranceDays_;
    Date today_;
};

class FmmLsmMultiLegOptionEngine : public MultiLegOption::engine {
public:
    FmmLsmMultiLegOptionEngine(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                               const QuantLib::ext::shared_ptr<FmmGrid>& grid, const FmmLsmEngineConfig& config);
    void calculate() const override;

private:
    QuantLib::ext::shared_ptr<ForwardMarketModel> model_;
    QuantLib::ext::shared_ptr<FmmGrid> grid_;
    FmmLsmEngineConfig config_;
};

//! writes the LSM statistics into an additional-results map under fmm* keys
void fmmWriteLsmResults(std::map<std::string, QuantLib::ext::any>& out, const FmmLsmResult& res,
                        const FmmDualBoundResult* dual, const FmmLsmEngineConfig& config, const FmmGrid& grid);

} // namespace QuantExt

#endif
