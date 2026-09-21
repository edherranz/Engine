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

#include <ql/any.hpp>
#include <ql/optional.hpp>
#include <ql/patterns/observable.hpp>

#include <functional>
#include <map>

namespace QuantExt {

//! exercise-policy treatment across revaluations (A6): retrain the LSM policy on every valuation
//! (retrained-policy risk) or keep the policy trained by the first valuation and only revalue
//! it on the same valuation paths (frozen-policy risk, paired common random numbers)
enum class FmmPolicyMode { Retrain, Frozen };

//! LSM engine configuration shared by the FMM trade engines
struct FmmLsmEngineConfig {
    FmmLsmConfig lsm;
    FmmPolicyMode policyMode = FmmPolicyMode::Retrain;
    bool dualBound = false; //!< also run the Andersen-Broadie dual bound (nested simulation)
    Size dualOuterPaths = 512, dualInnerPaths = 64;
    BigNatural dualSeed = 20260920;
    Natural gridToleranceDays = 3;
    //! observables the engine registers with (typically the calibrated IrModel adapter, so a
    //! recalibration invalidates cached instrument values)
    std::vector<QuantLib::ext::shared_ptr<Observable>> observables;
    //! calibration record merged into the additional results at pricing time (e.g. the
    //! builder's FmmCalibrationInfo::additionalResults)
    std::function<std::map<std::string, QuantLib::ext::any>()> calibrationResults;
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

//! ORE / QuantLib legs with signs (+1 receive, -1 pay from the holder's view) onto a vanilla
//! Cancel-style instrument without rights (a swap): the flow set the LSM helpers and the exposure
//! engine value exactly
FmmCallableInstrument fmmMapLegs(const std::vector<Leg>& legs, const std::vector<Real>& signs, const FmmGrid& grid,
                                 const Natural toleranceDays, const Date& today);

//! MultiLegOption arguments onto an Enter-style instrument from the option holder's view with one
//! right per future exercise (notice) date; ORE passes the notice dates as exercise dates and the
//! settlement (swap entry) dates alongside. Without an exercise object, or with no right left, the
//! returned instrument has no rights. usedSettle receives the settlement dates of the rights.
FmmCallableInstrument fmmMapMultiLegOption(const MultiLegOption::arguments& args, const FmmGrid& grid,
                                           const Natural toleranceDays, const Date& today,
                                           std::vector<Date>& usedSettle);

class FmmLsmMultiLegOptionEngine : public MultiLegOption::engine {
public:
    FmmLsmMultiLegOptionEngine(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                               const QuantLib::ext::shared_ptr<FmmGrid>& grid, const FmmLsmEngineConfig& config);
    void calculate() const override;

private:
    QuantLib::ext::shared_ptr<ForwardMarketModel> model_;
    QuantLib::ext::shared_ptr<FmmGrid> grid_;
    FmmLsmEngineConfig config_;
    mutable QuantLib::ext::optional<FmmLsmPolicy> frozenPolicy_;
};

//! writes the LSM statistics into an additional-results map under fmm* keys
void fmmWriteLsmResults(std::map<std::string, QuantLib::ext::any>& out, const FmmLsmResult& res,
                        const FmmDualBoundResult* dual, const FmmLsmEngineConfig& config, const FmmGrid& grid);

} // namespace QuantExt

#endif
