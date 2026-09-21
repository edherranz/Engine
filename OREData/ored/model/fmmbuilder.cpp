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

#include <ored/model/fmmbuilder.hpp>
#include <limits>

#include <ql/time/daycounters/actual365fixed.hpp>
#include <ored/model/structuredmodelerror.hpp>
#include <ored/model/structuredmodelwarning.hpp>
#include <ored/model/utilities.hpp>
#include <ored/utilities/log.hpp>
#include <ored/utilities/parsers.hpp>

#include <qle/models/fmmcalibration.hpp>
#include <qle/pricingengines/fmmanalyticswaptionengine.hpp>

#include <ql/math/solvers1d/brent.hpp>
#include <ql/models/shortrate/calibrationhelpers/swaptionhelper.hpp>
#include <ql/pricingengines/blackformula.hpp>

#include <algorithm>
#include <set>
#include <sstream>

using namespace QuantLib;
using namespace QuantExt;

namespace ore {
namespace data {

std::map<std::string, QuantLib::ext::any> FmmCalibrationInfo::additionalResults() const {
    std::map<std::string, QuantLib::ext::any> m;
    if (!valid)
        return m;
    // model-neutral keys (the LGM/HW builders emit the same information under lgm* names; the
    // FMM emits both a neutral and an fmm* spelling, never lgm* aliases)
    m["modelType"] = std::string("FMM");
    m["calibrationError"] = rmse;
    m["calibrationBasketExpiryTimes"] = expiryTimes;
    m["calibrationBasketStrikes"] = strikes;
    m["calibrationBasketMarketVols"] = marketVols;
    m["calibrationBasketModelVols"] = modelVols;
    m["fmmCalibrationError"] = rmse;
    m["fmmCalibrationApproximationMethod"] = approximationMethod;
    m["fmmCalibrationExpiryTimes"] = expiryTimes;
    m["fmmCalibrationStrikes"] = strikes;
    m["fmmCalibrationMarketVols"] = marketVols;
    m["fmmCalibrationModelVols"] = modelVols;
    m["fmmCalibrationEffectiveShifts"] = effectiveShifts;
    if (!mcVols.empty()) {
        m["fmmCalibrationMcVols"] = mcVols;
        m["fmmCalibrationMcVolStdErrors"] = mcVolStdErrors;
        m["fmmCalibrationMcIterations"] = static_cast<Real>(mcIterations);
        m["fmmCalibrationMcConverged"] = mcConverged;
    }
    if (!capletTimes.empty()) {
        m["fmmCalibrationCapletTimes"] = capletTimes;
        m["fmmCalibrationCapletStrikes"] = capletStrikes;
        m["fmmCalibrationCapletMarketVols"] = capletMarketVols;
        m["fmmCalibrationCapletModelVols"] = capletModelVols;
        m["fmmJointIterations"] = static_cast<Real>(jointIterations);
        m["fmmJointConverged"] = jointConverged;
        m["fmmJointStationary"] = jointStationary;
        m["fmmJointLastParameterChange"] = jointLastParameterChange;
        m["fmmCapletMaxResidualBp"] = capletMaxResidualBp;
    }
    m["fmmGridPeriods"] = static_cast<Real>(gridPeriods);
    m["fmmGridMaxDateMismatchDays"] = static_cast<Real>(gridMaxDateMismatchDays);
    std::vector<Real> t, len, k, fwd, ann, vega, sd;
    for (const auto& s : swaptionData) {
        t.push_back(s.timeToExpiry);
        len.push_back(s.swapLength);
        k.push_back(s.strike);
        fwd.push_back(s.atmForward);
        ann.push_back(s.annuity);
        vega.push_back(s.vega);
        sd.push_back(s.stdDev);
    }
    if (!swaptionData.empty()) {
        m["calibrationBasketSwapLengths"] = len;
        m["calibrationBasketAtmForwards"] = fwd;
        m["calibrationBasketAnnuities"] = ann;
        m["calibrationBasketVegas"] = vega;
        m["calibrationBasketStdDevs"] = sd;
    }
    return m;
}

namespace {
Real fmmGridTenorYears(const Period& g) {
    switch (g.units()) {
    case Days:
        return g.length() / 365.0;
    case Weeks:
        return 7.0 * g.length() / 365.0;
    case Months:
        return g.length() / 12.0;
    case Years:
        return static_cast<Real>(g.length());
    default:
        QL_FAIL("FmmBuilder: unsupported grid tenor units");
    }
}
// per-period scale of the level volatility: Lambda_j = tau_j / tau_grid * level, so that every
// period carries the same forward-rate volatility (HJM consistency); regular grids are unchanged
std::vector<Real> fmmLevelScale(const QuantExt::FmmGrid& grid, const Period& gridTenor) {
    const Real tauGrid = fmmGridTenorYears(gridTenor);
    std::vector<Real> scale(grid.numberOfRates());
    for (Size j = 1; j <= grid.numberOfRates(); ++j)
        scale[j - 1] = (grid.times()[j] - grid.times()[j - 1]) / tauGrid;
    return scale;
}
} // namespace

FmmBuilder::FmmBuilder(const QuantLib::ext::shared_ptr<ore::data::Market>& market,
                       const QuantLib::ext::shared_ptr<FmmData>& data, const std::string& configuration,
                       const Real bootstrapTolerance, const bool continueOnError,
                       const std::string& referenceCalibrationGrid, const bool setCalibrationInfo,
                       const std::string& id, BlackCalibrationHelper::CalibrationErrorType calibrationErrorType,
                       const bool allowChangingFallbacksUnderScenarios, const bool allowModelFallbacks,
                       const bool dontCalibrate)
    : IrModelBuilder(market, data, data->optionExpiries(), data->optionTerms(), data->optionStrikes(), configuration,
                     bootstrapTolerance, continueOnError, referenceCalibrationGrid, calibrationErrorType,
                     allowChangingFallbacksUnderScenarios, allowModelFallbacks,
                     data->calibrateVolatility() && data->calibrationType() != CalibrationType::None, dontCalibrate,
                     "FMM", id),
      setCalibrationInfo_(setCalibrationInfo) {
    data->validate();
    // the joint calibration also depends on the optionlet surface: observe it so that cap/floor
    // volatility scenarios trigger a recalibration (A6)
    if (data->capFloorBasket() == "ATM") {
        auto ovs = market->capFloorVol(data->qualifier(), configuration);
        if (!ovs.empty())
            marketObserver_->addObservable(*ovs);
    }
}

QuantLib::ext::shared_ptr<ForwardMarketModel> FmmBuilder::modelAsFmm() const {
    calculate();
    return fmm_;
}

QuantLib::ext::shared_ptr<FmmGrid> FmmBuilder::grid() const {
    calculate();
    return grid_;
}

QuantLib::ext::shared_ptr<FmmIrModel> FmmBuilder::irModel() const {
    calculate();
    return irModel_;
}

void FmmBuilder::initParametrization() const {
    if (parametrizationInitializedOnAnchorDate_ == referenceDate_)
        return;
    auto data = QuantLib::ext::dynamic_pointer_cast<FmmData>(data_);
    QL_REQUIRE(data, "FmmBuilder: data is not FmmData");

    // grid: the trade's contractual dates (when supplied) plus the calibration instruments' schedule
    // dates that are not within the mapping tolerance of a grid date already (helpers expiring on
    // notice dates carry schedules anchored on those dates, a month off the trade's lattice)
    std::set<Date> dates(data->gridDates().begin(), data->gridDates().end());
    std::set<Date> helperDates;
    for (const auto& h : swaptionBasket_) {
        auto sh = QuantLib::ext::dynamic_pointer_cast<SwaptionHelper>(h);
        QL_REQUIRE(sh, "FmmBuilder: calibration basket must consist of swaption helpers");
        const auto& swap = sh->underlying();
        helperDates.insert(swap->fixedSchedule().dates().begin(), swap->fixedSchedule().dates().end());
        helperDates.insert(swap->floatingSchedule().dates().begin(), swap->floatingSchedule().dates().end());
    }
    if (dates.empty())
        dates = helperDates;
    QL_REQUIRE(!dates.empty(), "FmmBuilder: no grid dates (neither trade dates nor a calibration basket) for "
                                   << data->qualifier());
    auto makeGrid = [&](const std::set<Date>& d) {
        return QuantLib::ext::make_shared<FmmGrid>(referenceDate_, d, data->grid(), Actual365Fixed(), 2,
                                                   data->noticePeriod(), data->noticeCalendar(),
                                                   data->noticeConvention());
    };
    grid_ = makeGrid(dates);
    Size added = 0;
    for (const Date& d : helperDates) {
        if (d <= referenceDate_)
            continue;
        auto it = std::lower_bound(grid_->dates().begin(), grid_->dates().end(), d);
        Integer best = std::numeric_limits<Integer>::max();
        if (it != grid_->dates().end())
            best = std::min(best, static_cast<Integer>(*it - d));
        if (it != grid_->dates().begin())
            best = std::min(best, static_cast<Integer>(d - *std::prev(it)));
        if (best > static_cast<Integer>(data->gridToleranceDays())) {
            dates.insert(d);
            ++added;
        }
    }
    if (added > 0) {
        DLOG("FMM grid for " << data->qualifier() << ": " << added
                             << " calibration instrument dates added beyond the mapping tolerance");
        grid_ = makeGrid(dates);
    }
    const Size M = grid_->numberOfRates();
    DLOG("FMM grid for " << data->qualifier() << ": " << M << " periods to " << grid_->dates().back());

    // volatility time grid: the swaption expiries (all but the last) for a bootstrap, else as given
    Array volTimes(data->volTimes().begin(), data->volTimes().end());
    std::vector<Real> levels = data->volValues();
    if (data->calibrateVolatility() && data->calibrationType() == CalibrationType::Bootstrap) {
        QL_REQUIRE(!swaptionExpiries_.empty() || data->capFloorBasket() != "None",
                   "FmmBuilder: no calibrating swaption provided.");
        if (swaptionExpiries_.empty()) {
            volTimes = Array(); // caplet-only bootstrap: constant a(t)
            levels = {data->volValues().front()};
        } else {
            volTimes = Array(swaptionExpiries_.begin(), std::next(swaptionExpiries_.end(), -1));
            levels.assign(volTimes.size() + 1, data->volValues().front());
        }
    } else if (data->calibrateVolatility() && data->calibrationType() == CalibrationType::BestFit &&
               data->bestFitTimeDependence()) {
        // joint / grid best fit: one a(t) segment per distinct swaption expiry (all but the last), the
        // duplicates of an expiry x term grid collapsed; a single expiry means a constant a(t)
        QL_REQUIRE(!swaptionExpiries_.empty(), "FmmBuilder: BestFitTimeDependence needs calibration swaptions.");
        std::vector<Time> distinct(swaptionExpiries_.begin(), swaptionExpiries_.end());
        std::sort(distinct.begin(), distinct.end());
        distinct.erase(std::unique(distinct.begin(), distinct.end(),
                                   [](const Time a, const Time b) { return std::fabs(a - b) < 1.0 / 365.0; }),
                       distinct.end());
        volTimes = distinct.size() <= 1 ? Array() : Array(distinct.begin(), std::next(distinct.end(), -1));
        levels.assign(volTimes.size() + 1, data->volValues().front());
    } else if (data->volParamType() == ParamType::Constant) {
        volTimes = Array();
        levels = {data->volValues().front()};
    } else if (levels.size() == 1) {
        levels.assign(volTimes.size() + 1, levels.front());
    }
    QL_REQUIRE(levels.size() == volTimes.size() + 1, "FmmBuilder: volatility values (" << levels.size()
                                                                                       << ") vs time grid ("
                                                                                       << volTimes.size() << ")");
    // the level volatility is quoted for a full grid-tenor period; a period of length tau carries
    // tau / tauGrid of it, so that every period has the same forward-rate volatility (HJM
    // consistency: Lambda_j = tau_j sigma_f; a flat Lambda over unequal periods would give a
    // one-month period three times the volatility of a quarterly one). Regular grids are unchanged.
    levelScale_ = fmmLevelScale(*grid_, data->grid());
    std::vector<Array> volLevels(M, Array(levels.begin(), levels.end()));
    Array shifts(M);
    for (Size j = 1; j <= M; ++j) {
        const Real tau = grid_->times()[j] - grid_->times()[j - 1];
        shifts[j - 1] = data->shift() == Null<Real>() ? 1.0 / tau : data->shift();
        volLevels[j - 1] *= levelScale_[j - 1];
    }

    fmmParametrization_ = QuantLib::ext::make_shared<FmmParametrization>(
        parseCurrency(currency_), modelDiscountCurve_, grid_->times(), shifts, volTimes, volLevels,
        data->volatilityType(), data->rhoInf(), data->beta(), data->factors());
    parametrization_ = fmmParametrization_;
    fmm_ = QuantLib::ext::make_shared<ForwardMarketModel>(fmmParametrization_);
    irModel_ = QuantLib::ext::make_shared<FmmIrModel>(fmm_);
    model_.linkTo(irModel_);
    params_[referenceDate_] = model_->params();
    parametrizationInitializedOnAnchorDate_ = referenceDate_;
}

void FmmBuilder::calibrate() const {
    auto data = QuantLib::ext::dynamic_pointer_cast<FmmData>(data_);
    const auto method = data->approximationMethod();
    auto& p = *fmmParametrization_;
    info_ = FmmCalibrationInfo();
    info_.approximationMethod = to_string(method);
    error_ = QL_MAX_REAL;
    const std::string errorTemplate =
        std::string("Failed to calibrate FMM Model. ") +
        (continueOnError_ ? std::string("Calculation will proceed.") : std::string("Calculation will be aborted."));
    try {
        // targets: each helper at its own strike, in normal-vol terms (a normal quote is taken as
        // is; a lognormal / shifted-lognormal quote is converted through the helper's market
        // price with the model's grid annuity)
        std::vector<FmmSwaptionVolTarget> targets;
        for (const auto& h : swaptionBasket_) {
            auto sh = QuantLib::ext::dynamic_pointer_cast<SwaptionHelper>(h);
            QL_REQUIRE(sh, "FmmBuilder: calibration basket must consist of swaption helpers");
            Swaption::arguments args;
            sh->swaption()->setupArguments(&args);
            Real K = 0.0, nominal = 1.0;
            Option::Type type = Option::Call;
            FmmSwaptionVolTarget t;
            t.swap = fmmSwapSpecFromSwaption(*grid_, args, data->gridToleranceDays(), K, type, nominal);
            t.strike = K;
            if (sh->volatilityType() == Normal) {
                t.normalVol = sh->volatility()->value();
            } else {
                const Real A0 = nominal * fmmAnnuity(p, t.swap), S0 = fmmForwardSwapRate(p, t.swap);
                t.normalVol = bachelierBlackFormulaImpliedVol(type, K, S0, p.rateTime(t.swap.a), sh->marketValue() / A0);
            }
            std::ostringstream lbl;
            lbl << "swaption_" << p.rateTime(t.swap.a) << "y_x_" << p.rateTime(t.swap.b) - p.rateTime(t.swap.a) << "y";
            t.label = lbl.str();
            targets.push_back(t);
        }
        std::sort(targets.begin(), targets.end(),
                  [](const FmmSwaptionVolTarget& a, const FmmSwaptionVolTarget& b) { return a.swap.a < b.swap.a; });

        FmmSeparableVols v;
        v.segmentTimes = p.parameterTimes(0);
        v.a.assign(v.segmentTimes.size() + 1, 1.0);
        QL_REQUIRE(levelScale_.size() == p.numberOfRates(), "FmmBuilder: level scale not initialised");
        auto scaledLevels = [&](const Real base) {
            std::vector<Real> L(p.numberOfRates());
            for (Size j = 0; j < L.size(); ++j)
                L[j] = base * levelScale_[j];
            return L;
        };
        v.levels = scaledLevels(data->volValues().front());
        // cap/floor basket: one backward-looking caplet per grid period up to the horizon, struck at
        // the period's forward, target vol from the market's stripped optionlet surface (looked up at
        // the period end, as ORE's overnight caplet pricer does)
        std::vector<FmmCapletVolTarget> caplets;
        if (data->capFloorBasket() == "ATM") {
            const auto ovs = market_->capFloorVol(data->qualifier(), configuration_);
            QL_REQUIRE(!ovs.empty(), "FmmBuilder: no cap/floor volatility surface for " << data->qualifier());
            QL_REQUIRE(ovs->volatilityType() == Normal,
                       "FmmBuilder: the cap/floor basket needs a normal optionlet surface for " << data->qualifier());
            const Date horizon = data->capFloorHorizon().length() > 0 ? referenceDate_ + data->capFloorHorizon()
                                                                        : grid_->dates().back();
            for (Size j = 1; j <= p.numberOfRates(); ++j) {
                if (grid_->dates()[j] > horizon || p.rateTime(j) <= 0.0)
                    continue;
                FmmCapletVolTarget c;
                c.bucket = j;
                c.backwardLooking = true;
                const Real F0 = (p.termStructure()->discount(p.rateTime(j - 1)) /
                                     p.termStructure()->discount(p.rateTime(j)) -
                                 1.0) /
                                p.tau(j);
                c.normalVol = ovs->volatility(grid_->dates()[j], F0, true);
                caplets.push_back(c);
                info_.capletTimes.push_back(p.rateTime(j));
                info_.capletStrikes.push_back(F0);
                info_.capletMarketVols.push_back(c.normalVol);
            }
            QL_REQUIRE(!caplets.empty(), "FmmBuilder: cap/floor basket is empty (horizon " << horizon << ")");
        }
        if (data->calibrationType() == CalibrationType::Bootstrap) {
            if (!caplets.empty() && targets.empty()) {
                fmmCapletLevelBootstrap(p, v, caplets);
                info_.jointIterations = 1;
                info_.jointConverged = true;
            } else if (!caplets.empty()) {
                QL_REQUIRE(targets.size() == v.a.size(),
                           "FmmBuilder: bootstrap needs one volatility segment per calibration swaption, got "
                               << v.a.size() << " segments for " << targets.size() << " active swaptions");
                const auto rep = fmmJointBootstrap(p, v, caplets, targets, data->jointMaxIterations(),
                                                   data->jointToleranceBp(), method);
                info_.jointIterations = rep.iterations;
                info_.jointConverged = rep.converged;
                info_.jointStationary = rep.stationary;
                info_.jointLastParameterChange = rep.lastParameterChange;
            } else {
                QL_REQUIRE(targets.size() == v.a.size(),
                           "FmmBuilder: bootstrap needs one volatility segment per calibration swaption, got "
                               << v.a.size() << " segments for " << targets.size() << " active swaptions");
                fmmSwaptionTimeDependenceBootstrap(p, v, targets, method);
            }
        } else if (data->calibrationType() == CalibrationType::BestFit && data->bestFitTimeDependence()) {
            // least-squares fit of the a(t) segments to the whole swaption basket (any expiry x term
            // set), jointly with the caplet levels when a cap/floor basket is configured
            QL_REQUIRE(!targets.empty(), "FmmBuilder: BestFitTimeDependence needs active calibration swaptions");
            if (!caplets.empty()) {
                const auto rep = fmmJointBestFit(p, v, caplets, targets, data->jointMaxIterations(),
                                                 data->jointToleranceBp(), method);
                info_.jointIterations = rep.iterations;
                info_.jointConverged = rep.converged;
                info_.jointStationary = rep.stationary;
                info_.jointLastParameterChange = rep.lastParameterChange;
            } else {
                fmmSwaptionTimeDependenceBestFit(p, v, targets, method);
            }
        } else if (data->calibrationType() == CalibrationType::BestFit) {
            // one common level, mean vol residual zero across the basket
            Brent solver;
            const Real level = solver.solve(
                [&](const Real l) {
                    v.levels = scaledLevels(l);
                    v.apply(p);
                    Real sum = 0.0;
                    for (const auto& t : targets)
                        sum += fmmSwaptionApprox(p, t.swap, fmmTargetStrike(p, t), Option::Call, method).normalVol -
                               t.normalVol;
                    return sum / static_cast<Real>(targets.size());
                },
                1e-12, data->volValues().front(), 1e-8, 50.0);
            v.levels = scaledLevels(level);
            v.apply(p);
        } else {
            QL_FAIL("FmmBuilder: unsupported calibration type " << data->calibrationType());
        }
        if (!caplets.empty()) {
            Real worst = 0.0;
            for (const auto& c : caplets) {
                const Real modelVol =
                    fmmCapletNormalVol(p, c.bucket, info_.capletStrikes[info_.capletModelVols.size()], c.backwardLooking);
                info_.capletModelVols.push_back(modelVol);
                worst = std::max(worst, std::fabs(modelVol - c.normalVol) * 1e4);
            }
            info_.capletMaxResidualBp = worst;
        }
        if (data->mcCorrection().enabled) {
            const auto& mc = data->mcCorrection();
            const auto rep = fmmMcCorrectedSwaptionBootstrap(p, fmm_, v, targets, method, mc.pathsPerRep, mc.reps,
                                                             mc.seed, mc.maxIterations, mc.toleranceBp);
            info_.mcIterations = rep.iterations;
            info_.mcConverged = rep.converged;
            for (const auto& r : rep.rows) {
                info_.mcVols.push_back(r.mcVol);
                info_.mcVolStdErrors.push_back(r.mcVolSe);
            }
        }
        for (const auto& t : targets) {
            const auto r = fmmSwaptionApprox(p, t.swap, fmmTargetStrike(p, t), Option::Call, method);
            info_.expiryTimes.push_back(p.rateTime(t.swap.a));
            info_.strikes.push_back(fmmTargetStrike(p, t));
            info_.marketVols.push_back(t.normalVol);
            info_.modelVols.push_back(r.normalVol);
            info_.effectiveShifts.push_back(r.effectiveShift == Null<Real>() ? 0.0 : r.effectiveShift);
        }
        // the helpers see the calibrated model through their (re-attached) engines
        irModel_->notifyObservers();
        for (auto& h : swaptionBasket_)
            h->update();
        error_ = getCalibrationError(swaptionBasket_);
    } catch (const std::exception& e) {
        StructuredModelErrorMessage(errorTemplate, e.what(), id_).log();
    }
    info_.rmse = error_;
    info_.gridPeriods = grid_->numberOfRates();
    info_.gridMaxDateMismatchDays = grid_->maxMismatchDays();
    const bool withinTolerance = std::fabs(error_) < bootstrapTolerance_ ||
                                 (data->calibrationType() == CalibrationType::BestFit && error_ != QL_MAX_REAL);
    if (withinTolerance) {
        if (Log::instance().filter(ORE_DEBUG) || setCalibrationInfo_) {
            try {
                DLOGGERSTREAM("Basket details:");
                DLOGGERSTREAM(getBasketDetails(info_.swaptionData));
            } catch (const std::exception& e) {
                WLOG("An error occurred: " << e.what());
            }
            DLOGGERSTREAM("rmse = " << error_);
        }
        info_.valid = true;
    } else {
        const std::string exceptionMessage = "FMM (" + data->qualifier() + ") calibration target function value (" +
                                             std::to_string(error_) + ") exceeds notification threshold (" +
                                             std::to_string(bootstrapTolerance_) + ").";
        StructuredModelWarningMessage(errorTemplate, exceptionMessage, id_).log();
        try {
            WLOGGERSTREAM(getBasketDetails(info_.swaptionData));
        } catch (const std::exception& e) {
            WLOG("An error occurred: " << e.what());
        }
        info_.valid = true;
        if (!continueOnError_)
            QL_FAIL(exceptionMessage);
    }
}

QuantLib::ext::shared_ptr<PricingEngine> FmmBuilder::getPricingEngine() const {
    auto data = QuantLib::ext::dynamic_pointer_cast<FmmData>(data_);
    return QuantLib::ext::make_shared<FmmAnalyticSwaptionEngine>(fmmParametrization_, grid_, data->approximationMethod(),
                                                                 data->gridToleranceDays());
}

} // namespace data
} // namespace ore
