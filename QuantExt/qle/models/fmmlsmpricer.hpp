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

/*! \file fmmlsmpricer.hpp
    \brief Longstaff-Schwartz pricer for on-grid Bermudan-callable structures under the FMM:
           independent training and valuation path sets with the exercise policy frozen before
           valuation (lower-bound estimator), exercise statistics as first-class outputs, imported
           policies / imported decisions for exercise-transfer scenarios, and the Andersen-Broadie
           dual upper bound by nested simulation. Owner rulings 2026-09-19 baked in.
    \ingroup models
*/

#ifndef quantext_fmm_lsm_pricer_hpp
#define quantext_fmm_lsm_pricer_hpp

#include <qle/models/forwardmarketmodel.hpp>

#include <functional>

namespace QuantExt {
using namespace QuantLib;

//! on-grid callable structure, unit notional, flows signed from the POLICY OWNER's perspective
/*! Cancel style: the owner holds the underlying flows and may terminate them: on exercise of a
    right, flows with grid index > settleIdx die, flows up to and including settleIdx are kept,
    and feeFlow (signed, from the owner's perspective) is added at T_settle. Enter style: no flows
    unless exercised; on exercise the owner receives the flows with grid index > settleIdx plus
    feeFlow at T_settle (physical settlement). */
struct FmmCallableInstrument {
    enum class Style { Enter, Cancel };
    Style style = Style::Cancel;
    Size lastFlowIdx = 0;             //!< largest grid index carrying a flow
    //! deterministic additive discounting spread over the RFR curve (FMM_SPEC.md section 8),
    //! applied to every flow and fee of this instrument: risky DF = P(t,T) exp(-spread (T-t))
    Real issuerSpread = 0.0;
    std::vector<Real> fixedFlows;     //!< size M+1; amount paid at T_j (index 0 unused)
    std::vector<Real> floatWeights;   //!< size M+1; coefficient of tau_j R_j(T_j) paid at T_j
    //! compounded-RFR coupon over the grid periods (startIdx, endIdx], paid at T_pay (payIdx >=
    //! endIdx): amount = weight * (prod_{k} (1 + tau_k R_k) - 1) + spreadAmount. Pathwise it equals
    //! weight * (B(T_end)/B(T_start) - 1) + spreadAmount by the bank-account identity, so the
    //! deflated sum telescopes exactly like grid-period floats (ore-fmm docs/A5_PLAN.md section 4)
    struct CompoundedFloat {
        Size startIdx = 0, endIdx = 0, payIdx = 0;
        Real weight = 0.0, spreadAmount = 0.0;
    };
    std::vector<CompoundedFloat> compoundedFloats;
    struct Right {
        Size noticeIdx = 0; //!< decision taken at T_notice (grid index)
        Size settleIdx = 0; //!< settlement/effective grid index, >= noticeIdx
        Real feeFlow = 0.0; //!< signed flow at T_settle when exercised (owner's perspective)
    };
    std::vector<Right> rights; //!< ascending notice indices
    //! consistency checks; requireRights = false admits a vanilla structure (no rights)
    void validate(const Size M, const bool requireRights = true) const;
};

//! t = 0 curve value of an instrument's flows (fixed, grid floats, compounded floats), with the
//! issuer-spread discounting; the underlying value of a Cancel-style structure
Real fmmUnderlyingCurveValue(const FmmCallableInstrument& instrument, const FmmParametrization& p);

struct FmmLsmConfig {
    Size trainingPaths = 16384;
    Size valuationPaths = 16384;
    BigNatural trainingSeed = 42;
    BigNatural valuationSeed = 424242;
    SequenceType trainingSequence = MersenneTwister;
    SequenceType valuationSequence = MersenneTwister;
    Size basisOrder = 2; //!< 1 or 2 (full second order incl. cross terms) on the 3 regressors
};

//! trained exercise policy: regression coefficients per right on the basis of the regressors
//! (coterminal-to-maturity par rate, front rate, deflated intrinsic exercise value). The
//! regression is fitted on, and the policy only ever exercises in, the in-the-money domain
//! (intrinsic value > 0) for both instrument styles.
struct FmmLsmPolicy {
    std::vector<Array> coefficients; //!< per right; empty entry = never exercise at that right
    Size basisOrder = 2;
};

//! Andersen-Broadie (2004) primal-dual result: the frozen-policy lower bound, the dual upper
//! bound from the policy-induced martingale estimated by nested simulation, and the gap
struct FmmDualBoundResult {
    Real lowerBound = 0.0, lowerBoundSe = 0.0; //!< policy value on the outer paths (small sample)
    //! duality gap E[max_r (h_r - M_r)] - V_0(policy) with the regret-only standard error: the
    //! outer-sample policy value enters every dual term as an exact constant (the term at the
    //! policy's own exercise date telescopes to it), so the gap is estimated far more precisely
    //! than the raw upper bound. Tight upper bound = a precise lower bound (e.g. the LSM
    //! valuation on 32k paths) + gap, standard errors combined in quadrature.
    Real gap = 0.0, gapSe = 0.0;
    //! lowerBound + gap on the outer sample; its s.e. (in quadrature) is dominated by the small-
    //! sample noise of lowerBound - prefer lower(precise) + gap
    Real upperBound = 0.0, upperBoundSe = 0.0;
    //! the same with the underlying-flow control variate (Cancel style; equal to the raw values
    //! for Enter): the realized deflated underlying total, whose mean is known exactly from the
    //! curve, is subtracted path by path (the gap is unchanged by it)
    Real lowerBoundCv = 0.0, lowerBoundCvSe = 0.0;
    Real upperBoundCv = 0.0, upperBoundCvSe = 0.0;
    Size outerPaths = 0, innerPaths = 0;
    Real runtimeSeconds = 0.0;
};

struct FmmLsmResult {
    Real lowerBound = 0.0;    //!< frozen-policy value on the independent valuation paths
    Real lowerBoundSe = 0.0;
    //! control-variate value (Cancel style): lowerBound - (underlyingValueMc - underlyingValue),
    //! i.e. the realized underlying total is replaced by its exact curve value path by path; the
    //! remaining noise is that of the exercise component only. Equal to lowerBound for Enter.
    Real lowerBoundCv = 0.0, lowerBoundCvSe = 0.0;
    Real underlyingValueMc = 0.0, underlyingValueMcSe = 0.0; //!< MC mean of the underlying total (Cancel)
    Real trainingValue = 0.0; //!< in-sample estimate (diagnostic only; upward-biased)
    std::vector<Real> exerciseProbability;   //!< per right, unconditional, valuation paths
    std::vector<Real> exerciseProbabilitySe;
    Real noExerciseProbability = 0.0;
    Real expectedExerciseTime = 0.0; //!< mean notice time conditional on exercise (0 if never)
    Real underlyingValue = 0.0;      //!< t=0 curve value of the underlying flows (Cancel) resp. 0
    Real runtimeSeconds = 0.0;
};

class FmmLsmPricer {
public:
    FmmLsmPricer(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model, const FmmCallableInstrument& instrument,
                 const FmmLsmConfig& config = FmmLsmConfig());

    //! train on the training paths, then value with the frozen policy on the valuation paths
    FmmLsmResult calculate();
    //! the policy trained by the last calculate() call
    const FmmLsmPolicy& policy() const { return policy_; }
    //! value THIS instrument under an imported (externally trained) policy on independent paths,
    //! i.e. the imported regression coefficients applied to THIS instrument's regressors; the
    //! imported policy's rights must align 1:1 with this instrument's rights. Meaningful when
    //! the regressors of the two instruments are comparable.
    FmmLsmResult valueWithPolicy(const FmmLsmPolicy& importedPolicy, const BigNatural seed) const;
    //! paired difference on COMMON valuation paths between own-policy and imported-policy values
    //! of this instrument: mean, standard error (the exercise-mismatch value with a paired CI)
    std::pair<Real, Real> pairedPolicyDifference(const FmmLsmPolicy& importedPolicy, const BigNatural seed) const;
    //! exercise-transfer scenario (configurable, not a claim about universal practice): THIS
    //! instrument is exercised whenever the policy OWNER's frozen policy exercises the owner's
    //! own instrument on the same path (the decisions are imported, evaluated on the owner's
    //! own regressors, e.g. the dealer cancels its hedge swap and the issuer follows by calling
    //! the note). Both pricers must share the model, the last flow index and the notice dates;
    //! the owner must have been trained (calculate()).
    FmmLsmResult valueWithImportedDecisions(const FmmLsmPricer& policyOwner, const BigNatural seed) const;
    //! paired difference on COMMON valuation paths: own-policy value minus imported-decision value
    //! of this instrument (mean, standard error)
    std::pair<Real, Real> pairedImportedDecisionDifference(const FmmLsmPricer& policyOwner,
                                                           const BigNatural seed) const;
    //! per-path deflated values of this instrument under the trained policy on the valuation
    //! paths of `seed` (two pricers on the same model, last flow index, sequence and seed see the
    //! SAME paths), with the realized deflated underlying totals (Cancel style; zeros for Enter)
    //! - for paired common-random-number comparisons and identities
    void pathValues(const BigNatural seed, std::vector<Real>& values, std::vector<Real>& underlyingTotals) const;
    //! Andersen-Broadie dual upper bound for the trained policy (A4 acceptance 4), nested
    //! simulation with innerPaths sub-paths at each right of each of outerPaths outer paths.
    //! Rights with a notice period (settleIdx > noticeIdx) use the adapted exercise payoff: the
    //! conditional expectation at the notice date of the deflated exercise cash flows, i.e. the
    //! flows between notice and settlement and the fee marked on the notice-date curve.
    FmmDualBoundResult dualBound(const Size outerPaths, const Size innerPaths, const BigNatural seed) const;

    //! Flow, marking and regressor helpers (pure functions of the instrument and the model
    //! state; public for the exposure engine, which reuses the pricer's conventions)
    //@{
    //! deterministic issuer-spread discount factor anchored at time 0
    Real spreadDf(const Time T) const { return std::exp(-instrument_.issuerSpread * T); }
    //! realized flow paid at T_j given the state at T_j (rates R_k, k <= j, are fixed)
    Real realizedFlow(const Size j, const ForwardMarketModel::State& state) const;
    //! deflated (by `bank`) value at the state's date T_atIdx of the flows paid at grid indices in
    //! (from, to], marked on the state's curve (compounded coupons partly accrued at T_atIdx use
    //! their fixed rates for the elapsed periods)
    Real markDeflated(const ForwardMarketModel::State& state, const Size atIdx, const Size from, const Size to,
                      const Real bank) const;
    //! deflated (by `bank`) value at the state's date of right r's fee paid at its settlement
    Real feeDeflated(const ForwardMarketModel::State& state, const Size r, const Real bank) const;
    //! regressors at a right given the state at its notice date
    Array regressorsAt(const ForwardMarketModel::State& state, const Size r, const Real bank) const;
    //! policy decision at right r for the given regressors
    bool exerciseDecision(const Array& x, const Size r, const FmmLsmPolicy& pol) const;
    //! basis functions of the three regressors (order 1: linear, order 2: full quadratic)
    Array basis(const Array& x) const;
    //@}

private:
    //! realized time-0-deflated policy value continuing from right `fromRight` along a fresh
    //! inner path started at `state` (grid index = notice of fromRight - 1 ... see .cpp),
    //! `pastFlows` = deflated flows already realized before the start index
    Real innerPolicyValue(const ForwardMarketModel::State& start, const Size startIdx, const Size fromRight,
                          const Real pastFlows, const FmmLsmPolicy& pol,
                          const std::vector<ForwardMarketModel::StepData>& steps,
                          const std::function<Real()>& normal) const;
    struct PathData {
        std::vector<Real> deflatedFlows;   // per grid index, flow_j / B(T_j)
        std::vector<Array> regressors;     // per right
        std::vector<Real> deflatedFee;     // per right, feeFlow / B(T_settle)
        std::vector<Real> switchedValue;   // per right, deflated value of flows with j > settleIdx
    };
    void simulate(const Size paths, const BigNatural seed, const SequenceType seq,
                  std::vector<PathData>& out) const;
    //! first right at which `pol` exercises along the path, -1 if never
    Integer decideRight(const PathData& p, const FmmLsmPolicy& pol) const;
    //! pathwise deflated value when exercising at right `exercisedRight` (-1 = never)
    Real pathValueAt(const PathData& p, const Integer exercisedRight) const;
    //! value and exercise statistics of the valuation paths given the exercised right per path
    FmmLsmResult summarize(const std::vector<PathData>& paths, const std::vector<Integer>& exercised) const;
    //! preconditions for common-path decision import from `policyOwner`
    void checkCommonPaths(const FmmLsmPricer& policyOwner, const char* where) const;

    QuantLib::ext::shared_ptr<ForwardMarketModel> model_;
    FmmCallableInstrument instrument_;
    FmmLsmConfig config_;
    FmmLsmPolicy policy_;
    Size M_;
    Real underlyingCurve_ = 0.0; //!< t=0 curve value of the underlying flows (Cancel), 0 for Enter
};

} // namespace QuantExt

#endif
