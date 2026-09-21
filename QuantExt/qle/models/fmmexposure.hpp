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

/*! \file qle/models/fmmexposure.hpp
    \brief FMM exposure engine (A7, route 1): conditional trade values of a jointly simulated
           netting set on exposure dates, for an external NPV cube consumed by ORE's post-processor
    \ingroup models
*/

#ifndef quantext_fmm_exposure_hpp
#define quantext_fmm_exposure_hpp

#include <qle/models/fmmlsmpricer.hpp>

#include <string>
#include <vector>

namespace QuantExt {
using namespace QuantLib;

//! one trade of the netting set: instrument on the common grid and its trained exercise policy
struct FmmExposureTrade {
    std::string id;
    FmmCallableInstrument instrument; //!< rights empty for a vanilla structure (swap, bond)
    FmmLsmPolicy policy;              //!< A4 policy (trained on the same model), ignored without rights
};

struct FmmExposureConfig {
    //! training paths: the conditional-value regressions per exposure date (the exercise policies
    //! are trained by the caller on the same model, see FmmLsmPricer); must differ from the
    //! exposure paths (seed)
    Size trainingPaths = 16384;
    BigNatural trainingSeed = 42;
    SequenceType trainingSequence = MersenneTwister;
    //! exposure paths: the cube
    Size exposurePaths = 10000;
    BigNatural exposureSeed = 424242;
    SequenceType exposureSequence = MersenneTwister;
    Size basisOrder = 2;            //!< 1 or 2 on the three regressors (as the A4 policies)
    Size minRegressionPaths = 100;  //!< below this many alive training paths the estimate is their mean
    //! flows paid exactly on an exposure date: excluded (paid) by default, ORE's valuation-date
    //! convention; true includes them in the value at that date
    bool includeFlowsOnExposureDate = false;
};

//! per trade and exposure date diagnostics
struct FmmExposureDiagnostics {
    Real aliveFraction = 0.0;    //!< exposure paths with a pending exercise decision (regression used)
    Real noticedFraction = 0.0;  //!< decided, settlement pending
    Real settledFraction = 0.0;  //!< exercised and settled (delivered / cancelled)
    Real exercisedFraction = 0.0; //!< cumulative exercise probability up to the date
    Size regressionPaths = 0;    //!< alive training paths the regression was fitted on
    Real regressionR2 = Null<Real>(); //!< of the regressed part (Cancel: the value over the exact underlying)
    Real regressionResidualSd = Null<Real>(); //!< in-sample, t-money
    //! out-of-sample bias on the alive exposure paths: mean of (estimate - realised future value)
    Real outOfSampleBias = Null<Real>(), outOfSampleBiasSe = Null<Real>();
    //! martingale identity in time-0 deflated units: mean over the exposure paths of
    //! V(t_d) / B(t_d) exp(-s t_d) + (deflated flows paid up to t_d under the policy); equals the
    //! time-zero value in expectation when the estimator is unbiased and no look-ahead leaks
    Real martingaleCheck = Null<Real>(), martingaleCheckSe = Null<Real>();
};

struct FmmExposureResult {
    std::vector<Time> times; //!< exposure times
    //! numeraire B(t_d) per exposure date and path, date-major: numeraire[d * paths + n]
    std::vector<Real> numeraire;
    //! deflated trade values V(t_d) / B(t_d) per trade, date-major as the numeraire
    std::vector<std::vector<Real>> deflatedValues;
    //! time-zero policy value per trade on the exposure paths (mean of the deflated policy flows)
    std::vector<Real> t0Values, t0ValuesSe;
    //! time-zero curve value of the underlying flows (the vanilla value; 0 for Enter style)
    std::vector<Real> t0CurveValues;
    std::vector<std::vector<FmmExposureDiagnostics>> diagnostics; //!< [trade][date]
    Size paths = 0;
    Real runtimeSeconds = 0.0;
};

/*! Conditional trade values on exposure dates from one FMM (A7 route 1).

    Every trade of the netting set is simulated on the same paths of the scenario model. Exercise
    decisions are taken on each path at the rights' notice dates with the trade's frozen A4 policy
    (notice-date state only). The value of a trade at an exposure date t_d on a path is
    - the exact curve value at t_d of the committed flows (from the path's state at t_d) once no
      exercise decision is pending: vanilla structures always, callable structures after the
      decision (Enter: delivered flows after settlement plus the fee at settlement; Cancel: the
      underlying flows up to settlement plus the fee; never exercised after the last notice:
      Enter 0, Cancel the remaining underlying flows) - notification does not erase the
      obligations up to settlement, delivered swaps stay until they mature;
    - otherwise (a decision is still pending) a conditional-value estimate: a least-squares
      regression, trained per exposure date on independent training paths, of the t_d-money value
      of the future flows under the policy on basis functions of the t_d state (coterminal par
      rate, front rate, exact curve value of the underlying flows after t_d), evaluated on the
      exposure path's own t_d state. For Cancel-style structures the exact value of the remaining
      underlying flows is taken out before the regression and added back (only the cancellation
      option is regressed), which keeps the estimate exact where no optionality is left. Flows on
      the exposure date follow the configured convention.
    Exposure times may lie between tenor dates (states at user times of the path generator) but not
    beyond the model grid. Values are in trade currency, deflated by the bank account B(t_d) (spot
    measure numeraire); the numeraire per date and path is returned for the aggregation data.
*/
class FmmExposureEngine {
public:
    FmmExposureEngine(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                      const std::vector<FmmExposureTrade>& trades, const std::vector<Time>& exposureTimes,
                      const FmmExposureConfig& config);
    FmmExposureResult calculate() const;

    /*! t-money value at the state's time of the instrument's flows with grid index in (from, to]
        that are still to be paid (T_j > t, or T_j >= t when includeFlowsAtT), marked on the
        state's curve: fixed flows and grid floats through the (extended) bonds, compounded floats
        with their fixed periods, the accruing period through the extended bond B(t) / B(T_{a-1});
        the issuer spread as a deterministic factor from t. Fees are not included. */
    static Real markValue(const ForwardMarketModel& model, const FmmCallableInstrument& inst,
                          const ForwardMarketModel::State& state, const Size from, const Size to,
                          const bool includeFlowsAtT);
    //! t-money value of right r's fee if its settlement is still to come
    static Real feeValue(const ForwardMarketModel& model, const FmmCallableInstrument& inst,
                         const ForwardMarketModel::State& state, const Size r, const bool includeFlowsAtT);

private:
    struct TradeData;
    struct PathTrade {
        std::vector<Real> policyFlows; //!< deflated (time 0) flows under the policy, per grid index
        Integer exercised = -1;
        Real total = 0.0;
    };
    void processPath(const FmmPathGenerator::Path& path, std::vector<PathTrade>& out) const;
    //! decision pending at time t on a path with the given exercised right
    bool decisionPending(const FmmCallableInstrument& inst, const Integer exercised, const Time t) const;
    //! exact t-money value of the committed flows and fees at the state's time
    Real committedValue(const TradeData& td, const PathTrade& pt, const ForwardMarketModel::State& st) const;
    Array regressors(const TradeData& td, const ForwardMarketModel::State& st) const;
    //! part of an alive trade's value that is exact on the state (Cancel: the remaining underlying flows)
    static Real exactPart(const FmmCallableInstrument& inst, const Array& x);
    /*! basis of the conditional-value regression: the front rate and the exact underlying value
        (the coterminal par rate is left out: it is nearly collinear with the underlying value, and
        the collinear pair makes the least-squares fit noisy in the tails, where PFE lives); order 1
        linear, order 2 with the squares and the cross term */
    Array exposureBasis(const Array& x) const;
    //! standardised least-squares fit of one exposure date (columns centred and scaled, tiny ridge)
    struct Fit {
        bool active = false;
        Real yMean = 0.0;
        Array beta, mean, scale;
        Real evaluate(const Array& b) const;
    };
    bool isDue(const Time T, const Time t) const;

    QuantLib::ext::shared_ptr<ForwardMarketModel> model_;
    std::vector<FmmExposureTrade> trades_;
    std::vector<Time> exposureTimes_;
    FmmExposureConfig config_;
    std::vector<QuantLib::ext::shared_ptr<FmmLsmPricer>> helpers_; //!< per trade: flow / regressor helpers
    std::vector<Time> userTimes_;
    std::vector<Size> gridPos_, expPos_; //!< positions of the grid times / exposure times in userTimes_
    Size lastFlow_ = 0;
};

} // namespace QuantExt

#endif
