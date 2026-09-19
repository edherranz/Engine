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
           valuation (lower-bound estimator), exercise statistics as first-class outputs, and an
           exported policy for imported-policy (exercise-transfer) scenarios. Owner rulings
           2026-09-19 baked in; Andersen-Broadie dual upper bound is a separate, explicitly
           outstanding A4 acceptance item.
    \ingroup models
*/

#ifndef quantext_fmm_lsm_pricer_hpp
#define quantext_fmm_lsm_pricer_hpp

#include <qle/models/forwardmarketmodel.hpp>

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
    std::vector<Real> fixedFlows;     //!< size M+1; amount paid at T_j (index 0 unused)
    std::vector<Real> floatWeights;   //!< size M+1; coefficient of tau_j R_j(T_j) paid at T_j
    struct Right {
        Size noticeIdx = 0; //!< decision taken at T_notice (grid index)
        Size settleIdx = 0; //!< settlement/effective grid index, >= noticeIdx
        Real feeFlow = 0.0; //!< signed flow at T_settle when exercised (owner's perspective)
    };
    std::vector<Right> rights; //!< ascending notice indices
    void validate(const Size M) const;
};

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
//! (coterminal-to-maturity par rate, front rate, deflated mark of the switched flows)
struct FmmLsmPolicy {
    std::vector<Array> coefficients; //!< per right; empty entry = never exercise at that right
    Size basisOrder = 2;
};

struct FmmLsmResult {
    Real lowerBound = 0.0;    //!< frozen-policy value on the independent valuation paths
    Real lowerBoundSe = 0.0;
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
    //! value THIS instrument under an imported (externally trained) policy on independent paths;
    //! the imported policy's rights must align 1:1 with this instrument's rights. Used for the
    //! exercise-transfer scenario (configurable, not a claim about universal practice).
    FmmLsmResult valueWithPolicy(const FmmLsmPolicy& importedPolicy, const BigNatural seed) const;
    //! paired difference on COMMON valuation paths between own-policy and imported-policy values
    //! of this instrument: mean, standard error (the exercise-mismatch value with a paired CI)
    std::pair<Real, Real> pairedPolicyDifference(const FmmLsmPolicy& importedPolicy, const BigNatural seed) const;

private:
    struct PathData {
        std::vector<Real> deflatedFlows;   // per grid index, flow_j / B(T_j)
        std::vector<Array> regressors;     // per right
        std::vector<Real> deflatedFee;     // per right, feeFlow / B(T_settle)
        std::vector<Real> switchedValue;   // per right, deflated value of flows with j > settleIdx
    };
    void simulate(const Size paths, const BigNatural seed, const SequenceType seq,
                  std::vector<PathData>& out) const;
    Array basis(const Array& x) const;
    //! pathwise value under a policy (forward pass); optionally records the exercised right
    Real pathValue(const PathData& p, const FmmLsmPolicy& pol, Integer* exercisedRight) const;

    QuantLib::ext::shared_ptr<ForwardMarketModel> model_;
    FmmCallableInstrument instrument_;
    FmmLsmConfig config_;
    FmmLsmPolicy policy_;
    Size M_;
};

} // namespace QuantExt

#endif
