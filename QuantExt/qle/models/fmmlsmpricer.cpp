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

#include <qle/models/fmmlsmpricer.hpp>

#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/matrixutilities/qrdecomposition.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
#include <ql/math/statistics/incrementalstatistics.hpp>

#include <chrono>

namespace QuantExt {

void FmmCallableInstrument::validate(const Size M) const {
    QL_REQUIRE(lastFlowIdx >= 1 && lastFlowIdx <= M, "FmmCallableInstrument: lastFlowIdx out of range");
    QL_REQUIRE(fixedFlows.size() == M + 1 && floatWeights.size() == M + 1,
               "FmmCallableInstrument: flow vectors must have size M+1");
    QL_REQUIRE(!rights.empty(), "FmmCallableInstrument: no exercise rights");
    for (Size r = 0; r < rights.size(); ++r) {
        QL_REQUIRE(rights[r].noticeIdx >= 1 && rights[r].noticeIdx < lastFlowIdx,
                   "FmmCallableInstrument: notice index " << rights[r].noticeIdx << " out of range");
        QL_REQUIRE(rights[r].settleIdx >= rights[r].noticeIdx && rights[r].settleIdx <= lastFlowIdx,
                   "FmmCallableInstrument: settle index " << rights[r].settleIdx << " invalid");
        QL_REQUIRE(r == 0 || rights[r].noticeIdx > rights[r - 1].noticeIdx,
                   "FmmCallableInstrument: notice indices must be strictly ascending");
    }
}

FmmLsmPricer::FmmLsmPricer(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                           const FmmCallableInstrument& instrument, const FmmLsmConfig& config)
    : model_(model), instrument_(instrument), config_(config), M_(model->parametrization()->numberOfRates()) {
    QL_REQUIRE(model_, "FmmLsmPricer: model is null");
    instrument_.validate(M_);
    QL_REQUIRE(config_.basisOrder == 1 || config_.basisOrder == 2, "FmmLsmPricer: basisOrder must be 1 or 2");
    if (instrument_.style == FmmCallableInstrument::Style::Cancel) {
        const auto& p = *model_->parametrization();
        for (Size j = 1; j <= instrument_.lastFlowIdx; ++j) {
            const Real pj = p.termStructure()->discount(p.rateTime(j));
            const Real pjm = p.termStructure()->discount(p.rateTime(j - 1));
            underlyingCurve_ += (instrument_.fixedFlows[j] * pj + instrument_.floatWeights[j] * (pjm - pj)) *
                                spreadDf(p.rateTime(j));
        }
    }
}

Array FmmLsmPricer::basis(const Array& x) const {
    if (config_.basisOrder == 1) {
        Array b(4);
        b[0] = 1.0;
        b[1] = x[0];
        b[2] = x[1];
        b[3] = x[2];
        return b;
    }
    Array b(10);
    b[0] = 1.0;
    b[1] = x[0];
    b[2] = x[1];
    b[3] = x[2];
    b[4] = x[0] * x[0];
    b[5] = x[1] * x[1];
    b[6] = x[2] * x[2];
    b[7] = x[0] * x[1];
    b[8] = x[0] * x[2];
    b[9] = x[1] * x[2];
    return b;
}

Real FmmLsmPricer::markDeflated(const ForwardMarketModel::State& st, const Size from, const Size to,
                                const Real bank) const {
    // time-0 deflated conditional expectation of the flows: E_t[flow_j / B(T_j)] = flow_j P(t,T_j) / B(t),
    // times the deterministic issuer-spread factor to T_j - the same units as the realized
    // deflated flows (a time-t-relative spread factor here would be inconsistent with them)
    const auto& p = *model_->parametrization();
    Real mark = 0.0;
    for (Size j = from + 1; j <= to; ++j) {
        const Real pj = model_->discountBond(st, p.rateTime(j));
        const Real pjm = model_->discountBond(st, p.rateTime(j - 1));
        mark += (instrument_.fixedFlows[j] * pj + instrument_.floatWeights[j] * (pjm - pj)) *
                spreadDf(p.rateTime(j));
    }
    return mark / bank;
}

Real FmmLsmPricer::feeDeflated(const ForwardMarketModel::State& st, const Size r, const Real bank) const {
    const auto& p = *model_->parametrization();
    const auto& rt = instrument_.rights[r];
    return rt.feeFlow * model_->discountBond(st, p.rateTime(rt.settleIdx)) * spreadDf(p.rateTime(rt.settleIdx)) /
           bank;
}

Array FmmLsmPricer::regressorsAt(const ForwardMarketModel::State& st, const Size r, const Real bank) const {
    const auto& p = *model_->parametrization();
    const auto& rt = instrument_.rights[r];
    const Size last = instrument_.lastFlowIdx;
    // par rate to maturity, front rate, deflated intrinsic exercise value (what exercising is
    // worth over holding, marked on the current curve: Enter switched flows + fee, Cancel fee -
    // switched flows)
    Real annuity = 0.0;
    const Real pEnd = model_->discountBond(st, p.rateTime(last));
    for (Size j = rt.noticeIdx + 1; j <= last; ++j)
        annuity += p.tau(j) * model_->discountBond(st, p.rateTime(j));
    const Real mark = markDeflated(st, rt.settleIdx, last, bank);
    const Real fee = feeDeflated(st, r, bank);
    Array x(3);
    x[0] = annuity > QL_EPSILON ? (1.0 - pEnd) / annuity : 0.0;
    x[1] = st.R[std::min<Size>(rt.noticeIdx + 1, M_) - 1];
    x[2] = instrument_.style == FmmCallableInstrument::Style::Enter ? mark + fee : fee - mark;
    return x;
}

bool FmmLsmPricer::exerciseDecision(const Array& x, const Size r, const FmmLsmPolicy& pol) const {
    if (pol.coefficients[r].empty())
        return false;
    if (x[2] <= 0.0)
        return false; // no intrinsic value: outside the trained (in-the-money) domain
    const Array b = basis(x);
    Real dHat = 0.0;
    for (Size k = 0; k < b.size(); ++k)
        dHat += b[k] * pol.coefficients[r][k];
    return dHat < 0.0;
}

void FmmLsmPricer::simulate(const Size paths, const BigNatural seed, const SequenceType seq,
                            std::vector<PathData>& out) const {
    const auto& p = *model_->parametrization();
    const Size last = instrument_.lastFlowIdx;
    std::vector<Time> userTimes;
    for (Size j = 1; j <= last; ++j)
        userTimes.push_back(p.rateTime(j));
    FmmPathGenerator gen(model_, userTimes, seq, seed);
    const Size nRights = instrument_.rights.size();
    out.assign(paths, PathData());
    for (Size n = 0; n < paths; ++n) {
        const auto path = gen.next();
        PathData& d = out[n];
        d.deflatedFlows.assign(last + 1, 0.0);
        Real total = 0.0;
        std::vector<Real> bank(last + 1, 1.0);
        for (Size j = 1; j <= last; ++j) {
            const auto& st = path.states[j - 1];
            bank[j] = model_->bankAccount(st);
            const Real flow = instrument_.fixedFlows[j] + instrument_.floatWeights[j] * p.tau(j) * st.R[j - 1];
            d.deflatedFlows[j] = flow / bank[j] * spreadDf(p.rateTime(j));
            total += d.deflatedFlows[j];
        }
        d.regressors.resize(nRights);
        d.deflatedFee.resize(nRights);
        d.switchedValue.resize(nRights);
        for (Size r = 0; r < nRights; ++r) {
            const auto& rt = instrument_.rights[r];
            const auto& st = path.states[rt.noticeIdx - 1]; // state at T_notice
            Real sw = 0.0;
            for (Size j = rt.settleIdx + 1; j <= last; ++j)
                sw += d.deflatedFlows[j];
            d.switchedValue[r] = sw;
            d.deflatedFee[r] = feeDeflated(st, r, bank[rt.noticeIdx]);
            d.regressors[r] = regressorsAt(st, r, bank[rt.noticeIdx]);
        }
        d.deflatedFlows[0] = total; // cache the never-exercise total in slot 0
    }
}

Integer FmmLsmPricer::decideRight(const PathData& d, const FmmLsmPolicy& pol) const {
    for (Size r = 0; r < instrument_.rights.size(); ++r)
        if (exerciseDecision(d.regressors[r], r, pol))
            return static_cast<Integer>(r);
    return -1;
}

Real FmmLsmPricer::pathValueAt(const PathData& d, const Integer exercisedRight) const {
    const bool enter = instrument_.style == FmmCallableInstrument::Style::Enter;
    const Real total = d.deflatedFlows[0];
    if (exercisedRight < 0)
        return enter ? 0.0 : total;
    const Size r = static_cast<Size>(exercisedRight);
    return enter ? d.switchedValue[r] + d.deflatedFee[r] : total - d.switchedValue[r] + d.deflatedFee[r];
}

FmmLsmResult FmmLsmPricer::calculate() {
    const auto start = std::chrono::steady_clock::now();
    const bool enter = instrument_.style == FmmCallableInstrument::Style::Enter;
    const Size nRights = instrument_.rights.size();

    // ---- training pass ----
    std::vector<PathData> train;
    simulate(config_.trainingPaths, config_.trainingSeed, config_.trainingSequence, train);
    policy_.basisOrder = config_.basisOrder;
    policy_.coefficients.assign(nRights, Array());
    std::vector<Real> V(train.size());
    for (Size n = 0; n < train.size(); ++n)
        V[n] = enter ? 0.0 : train[n].deflatedFlows[0];
    const Size nb = basis(Array(3, 0.0)).size();
    for (Size rr = nRights; rr-- > 0;) {
        // regress on the in-the-money paths only (positive intrinsic value), both styles
        std::vector<Size> idx;
        for (Size n = 0; n < train.size(); ++n)
            if (train[n].regressors[rr][2] > 0.0)
                idx.push_back(n);
        if (idx.size() < 2 * nb)
            continue; // leave as never-exercise
        Matrix A(idx.size(), nb);
        Array y(idx.size());
        for (Size q = 0; q < idx.size(); ++q) {
            const PathData& d = train[idx[q]];
            const Array b = basis(d.regressors[rr]);
            for (Size k = 0; k < nb; ++k)
                A[q][k] = b[k];
            y[q] = V[idx[q]] - pathValueAt(d, static_cast<Integer>(rr)); // continuation minus exercise
        }
        policy_.coefficients[rr] = qrSolve(A, y);
        for (const Size n : idx) {
            const PathData& d = train[n];
            if (exerciseDecision(d.regressors[rr], rr, policy_))
                V[n] = pathValueAt(d, static_cast<Integer>(rr));
        }
    }
    Real trainMean = 0.0;
    for (const Real v : V)
        trainMean += v;
    trainMean /= static_cast<Real>(V.size());

    // ---- independent valuation pass with the frozen policy ----
    FmmLsmResult res = valueWithPolicy(policy_, config_.valuationSeed);
    res.trainingValue = trainMean;
    res.runtimeSeconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() *
        1e-3;
    return res;
}

FmmLsmResult FmmLsmPricer::summarize(const std::vector<PathData>& val, const std::vector<Integer>& exercised) const {
    const auto& p = *model_->parametrization();
    const bool enter = instrument_.style == FmmCallableInstrument::Style::Enter;
    const Size nRights = instrument_.rights.size();
    FmmLsmResult res;
    IncrementalStatistics stats, cvStats, totalStats, neverStats;
    std::vector<IncrementalStatistics> exStats(nRights);
    Real timeSum = 0.0;
    Size timeCount = 0;
    for (Size n = 0; n < val.size(); ++n) {
        const Integer er = exercised[n];
        const Real v = pathValueAt(val[n], er);
        stats.add(v);
        // control variate: the realized underlying total has the known mean underlyingCurve_
        const Real total = enter ? 0.0 : val[n].deflatedFlows[0];
        totalStats.add(total);
        cvStats.add(v - (total - underlyingCurve_));
        for (Size r = 0; r < nRights; ++r)
            exStats[r].add(er == static_cast<Integer>(r) ? 1.0 : 0.0);
        neverStats.add(er < 0 ? 1.0 : 0.0);
        if (er >= 0) {
            timeSum += p.rateTime(instrument_.rights[er].noticeIdx);
            ++timeCount;
        }
    }
    res.lowerBound = stats.mean();
    res.lowerBoundSe = stats.errorEstimate();
    res.lowerBoundCv = cvStats.mean();
    res.lowerBoundCvSe = cvStats.errorEstimate();
    res.underlyingValue = underlyingCurve_;
    res.underlyingValueMc = totalStats.mean();
    res.underlyingValueMcSe = totalStats.errorEstimate();
    for (Size r = 0; r < nRights; ++r) {
        res.exerciseProbability.push_back(exStats[r].mean());
        res.exerciseProbabilitySe.push_back(exStats[r].errorEstimate());
    }
    res.noExerciseProbability = neverStats.mean();
    res.expectedExerciseTime = timeCount > 0 ? timeSum / static_cast<Real>(timeCount) : 0.0;
    return res;
}

FmmLsmResult FmmLsmPricer::valueWithPolicy(const FmmLsmPolicy& pol, const BigNatural seed) const {
    QL_REQUIRE(pol.coefficients.size() == instrument_.rights.size(),
               "FmmLsmPricer::valueWithPolicy: policy rights (" << pol.coefficients.size()
                                                                << ") do not align with instrument rights ("
                                                                << instrument_.rights.size() << ")");
    QL_REQUIRE(pol.basisOrder == config_.basisOrder, "FmmLsmPricer::valueWithPolicy: basis order mismatch");
    std::vector<PathData> val;
    simulate(config_.valuationPaths, seed, config_.valuationSequence, val);
    std::vector<Integer> ex(val.size());
    for (Size n = 0; n < val.size(); ++n)
        ex[n] = decideRight(val[n], pol);
    return summarize(val, ex);
}

void FmmLsmPricer::pathValues(const BigNatural seed, std::vector<Real>& values,
                              std::vector<Real>& underlyingTotals) const {
    QL_REQUIRE(policy_.coefficients.size() == instrument_.rights.size(),
               "FmmLsmPricer::pathValues: no trained policy - call calculate() first");
    const bool enter = instrument_.style == FmmCallableInstrument::Style::Enter;
    std::vector<PathData> val;
    simulate(config_.valuationPaths, seed, config_.valuationSequence, val);
    values.resize(val.size());
    underlyingTotals.resize(val.size());
    for (Size n = 0; n < val.size(); ++n) {
        values[n] = pathValueAt(val[n], decideRight(val[n], policy_));
        underlyingTotals[n] = enter ? 0.0 : val[n].deflatedFlows[0];
    }
}

std::pair<Real, Real> FmmLsmPricer::pairedPolicyDifference(const FmmLsmPolicy& importedPolicy,
                                                           const BigNatural seed) const {
    QL_REQUIRE(importedPolicy.coefficients.size() == instrument_.rights.size(),
               "FmmLsmPricer::pairedPolicyDifference: policy rights do not align");
    std::vector<PathData> val;
    simulate(config_.valuationPaths, seed, config_.valuationSequence, val);
    IncrementalStatistics diff;
    for (const auto& d : val)
        diff.add(pathValueAt(d, decideRight(d, policy_)) - pathValueAt(d, decideRight(d, importedPolicy)));
    return {diff.mean(), diff.errorEstimate()};
}

void FmmLsmPricer::checkCommonPaths(const FmmLsmPricer& owner, const char* where) const {
    QL_REQUIRE(owner.model_ == model_, where << ": the policy owner must be priced on the same model (common paths)");
    QL_REQUIRE(owner.instrument_.lastFlowIdx == instrument_.lastFlowIdx,
               where << ": the policy owner's last flow index must agree (common paths)");
    QL_REQUIRE(owner.instrument_.rights.size() == instrument_.rights.size(), where << ": rights do not align 1:1");
    for (Size r = 0; r < instrument_.rights.size(); ++r)
        QL_REQUIRE(owner.instrument_.rights[r].noticeIdx == instrument_.rights[r].noticeIdx,
                   where << ": notice dates of right " << r << " differ");
    QL_REQUIRE(owner.policy_.coefficients.size() == owner.instrument_.rights.size(),
               where << ": the policy owner has no trained policy - call calculate() on it first");
}

FmmLsmResult FmmLsmPricer::valueWithImportedDecisions(const FmmLsmPricer& owner, const BigNatural seed) const {
    checkCommonPaths(owner, "FmmLsmPricer::valueWithImportedDecisions");
    std::vector<PathData> mine, theirs;
    simulate(config_.valuationPaths, seed, config_.valuationSequence, mine);
    owner.simulate(config_.valuationPaths, seed, config_.valuationSequence, theirs);
    std::vector<Integer> ex(mine.size());
    for (Size n = 0; n < mine.size(); ++n)
        ex[n] = owner.decideRight(theirs[n], owner.policy_);
    return summarize(mine, ex);
}

std::pair<Real, Real> FmmLsmPricer::pairedImportedDecisionDifference(const FmmLsmPricer& owner,
                                                                     const BigNatural seed) const {
    checkCommonPaths(owner, "FmmLsmPricer::pairedImportedDecisionDifference");
    std::vector<PathData> mine, theirs;
    simulate(config_.valuationPaths, seed, config_.valuationSequence, mine);
    owner.simulate(config_.valuationPaths, seed, config_.valuationSequence, theirs);
    IncrementalStatistics diff;
    for (Size n = 0; n < mine.size(); ++n)
        diff.add(pathValueAt(mine[n], decideRight(mine[n], policy_)) -
                 pathValueAt(mine[n], owner.decideRight(theirs[n], owner.policy_)));
    return {diff.mean(), diff.errorEstimate()};
}

Real FmmLsmPricer::innerPolicyValue(const ForwardMarketModel::State& start, const Size startIdx,
                                    const Size fromRight, const Real pastFlows, const FmmLsmPolicy& pol,
                                    const std::vector<ForwardMarketModel::StepData>& steps,
                                    const std::function<Real()>& normal) const {
    const bool enter = instrument_.style == FmmCallableInstrument::Style::Enter;
    const auto& p = *model_->parametrization();
    const Size last = instrument_.lastFlowIdx;
    const Size nRights = instrument_.rights.size();
    ForwardMarketModel::State st = start;
    Real acc = pastFlows;      // Cancel: all realized deflated flows so far
    Real switched = 0.0;       // Enter: deflated flows received after exercise
    bool exercised = false;
    Size settleIdx = 0;
    Real fee = 0.0;
    Size r = fromRight;
    Array z(M_ + 1);
    for (Size j = startIdx + 1; j <= last; ++j) {
        for (Size q = 0; q <= M_; ++q)
            z[q] = normal();
        model_->evolve(st, steps[j - 1], z);
        const Real bank = model_->bankAccount(st);
        const Real flow = (instrument_.fixedFlows[j] + instrument_.floatWeights[j] * p.tau(j) * st.R[j - 1]) /
                          bank * spreadDf(p.rateTime(j));
        if (exercised) {
            if (j > settleIdx)
                switched += flow;
            else if (!enter)
                acc += flow;
            continue;
        }
        acc += flow;
        if (r < nRights && instrument_.rights[r].noticeIdx == j) {
            const Array x = regressorsAt(st, r, bank);
            if (exerciseDecision(x, r, pol)) {
                exercised = true;
                settleIdx = instrument_.rights[r].settleIdx;
                fee = feeDeflated(st, r, bank);
                if (!enter && settleIdx == j)
                    return acc + fee; // Cancel with immediate settlement: value is known
            }
            ++r;
        }
    }
    if (!exercised)
        return enter ? 0.0 : acc;
    return enter ? switched + fee : acc + fee;
}

FmmDualBoundResult FmmLsmPricer::dualBound(const Size outerPaths, const Size innerPaths,
                                           const BigNatural seed) const {
    const auto start = std::chrono::steady_clock::now();
    QL_REQUIRE(policy_.coefficients.size() == instrument_.rights.size(),
               "FmmLsmPricer::dualBound: no trained policy - call calculate() first");
    const bool enter = instrument_.style == FmmCallableInstrument::Style::Enter;
    const auto& p = *model_->parametrization();
    const Size last = instrument_.lastFlowIdx;
    const Size nRights = instrument_.rights.size();

    std::vector<ForwardMarketModel::StepData> steps;
    for (Size j = 1; j <= last; ++j)
        steps.push_back(model_->makeStep(p.rateTime(j - 1), p.rateTime(j)));

    // outer paths: states at every grid date, deflated flows, prefix sums
    std::vector<Time> userTimes;
    for (Size j = 1; j <= last; ++j)
        userTimes.push_back(p.rateTime(j));
    FmmPathGenerator gen(model_, userTimes, MersenneTwister, seed);
    struct Outer {
        std::vector<ForwardMarketModel::State> states;
        std::vector<Real> bank, prefix; // prefix[j] = sum of deflated flows up to j
        Real lower = 0.0;
        Real cv = 0.0; // realized underlying total minus its curve value (Cancel), 0 for Enter
    };
    std::vector<Outer> outer(outerPaths);
    IncrementalStatistics lowerStats, lowerCvStats;
    for (Size n = 0; n < outerPaths; ++n) {
        const auto path = gen.next();
        Outer& o = outer[n];
        o.states = path.states;
        o.bank.assign(last + 1, 1.0);
        o.prefix.assign(last + 1, 0.0);
        for (Size j = 1; j <= last; ++j) {
            o.bank[j] = model_->bankAccount(o.states[j - 1]);
            const Real flow =
                (instrument_.fixedFlows[j] + instrument_.floatWeights[j] * p.tau(j) * o.states[j - 1].R[j - 1]) /
                o.bank[j] * spreadDf(p.rateTime(j));
            o.prefix[j] = o.prefix[j - 1] + flow;
        }
        // frozen-policy value along the outer path (lower-bound sample)
        Real lower = enter ? 0.0 : o.prefix[last];
        for (Size r = 0; r < nRights; ++r) {
            const auto& rt = instrument_.rights[r];
            const auto& st = o.states[rt.noticeIdx - 1];
            if (exerciseDecision(regressorsAt(st, r, o.bank[rt.noticeIdx]), r, policy_)) {
                const Real fee = feeDeflated(st, r, o.bank[rt.noticeIdx]);
                lower = enter ? (o.prefix[last] - o.prefix[rt.settleIdx]) + fee : o.prefix[rt.settleIdx] + fee;
                break;
            }
        }
        o.lower = lower;
        o.cv = enter ? 0.0 : o.prefix[last] - underlyingCurve_;
        lowerStats.add(lower);
        lowerCvStats.add(lower - o.cv);
    }
    const Real lowerGlobal = lowerStats.mean();

    // dual pass: martingale from the policy value function, nested continuation estimates
    MersenneTwisterUniformRng rng(seed + 7919);
    InverseCumulativeNormal icn;
    auto normal = [&rng, &icn]() { return icn(rng.next().value); };
    IncrementalStatistics upperStats, upperCvStats;
    for (Size n = 0; n < outerPaths; ++n) {
        const Outer& o = outer[n];
        Real M = 0.0, maxTerm = -QL_MAX_REAL;
        Real prevExpectation = lowerGlobal; // E[L_0] estimate
        Real contHat = 0.0;
        for (Size r = 0; r < nRights; ++r) {
            const auto& rt = instrument_.rights[r];
            const auto& st = o.states[rt.noticeIdx - 1];
            const Real bank = o.bank[rt.noticeIdx];
            const Array x = regressorsAt(st, r, bank);
            // adapted exercise payoff: the conditional expectation at the notice date of the
            // deflated exercise cash flows - Enter: intrinsic value (switched flows + fee marked
            // on the notice-date curve); Cancel: realized flows to notice, plus the flows between
            // notice and settlement and the fee marked on the notice-date curve. Realized post-
            // notice flows would leak information and invalidate the dual bound.
            const Real h = enter ? x[2]
                                 : o.prefix[rt.noticeIdx] + markDeflated(st, rt.noticeIdx, rt.settleIdx, bank) +
                                       feeDeflated(st, r, bank);
            // continuation under the policy from the next right: nested inner simulation
            if (enter && r + 1 == nRights) {
                contHat = 0.0; // nothing left to enter
            } else {
                IncrementalStatistics inner;
                for (Size m = 0; m < innerPaths; ++m)
                    inner.add(innerPolicyValue(st, rt.noticeIdx, r + 1, o.prefix[rt.noticeIdx], policy_, steps,
                                               normal));
                contHat = inner.mean();
            }
            const bool ex = exerciseDecision(x, r, policy_);
            const Real L = ex ? h : contHat;
            M += L - prevExpectation;
            maxTerm = std::max(maxTerm, h - M);
            prevExpectation = contHat; // E_r[L_{r+1}] shares the same estimator
        }
        // terminal (never exercise): realized terminal value
        const Real hT = enter ? 0.0 : o.prefix[last];
        M += hT - prevExpectation;
        maxTerm = std::max(maxTerm, hT - M);
        upperStats.add(maxTerm);
        upperCvStats.add(maxTerm - o.cv);
    }

    // max_r (h_r - M_r) = V^_0 + regret on every outer path (V^_0 = lowerGlobal enters each term
    // as an exact constant; the term at the policy's own exercise date telescopes to it), so the
    // duality gap mean(maxTerm) - V^_0 carries only the regret noise, whereas the raw upper bound
    // inherits the full small-sample noise of V^_0
    FmmDualBoundResult res;
    res.lowerBound = lowerGlobal;
    res.lowerBoundSe = lowerStats.errorEstimate();
    res.gap = upperStats.mean() - lowerGlobal;
    res.gapSe = upperStats.errorEstimate();
    res.upperBound = upperStats.mean();
    res.upperBoundSe = std::sqrt(res.lowerBoundSe * res.lowerBoundSe + res.gapSe * res.gapSe);
    res.lowerBoundCv = lowerCvStats.mean();
    res.lowerBoundCvSe = lowerCvStats.errorEstimate();
    res.upperBoundCv = upperCvStats.mean();
    res.upperBoundCvSe = std::sqrt(res.lowerBoundCvSe * res.lowerBoundCvSe + res.gapSe * res.gapSe);
    res.outerPaths = outerPaths;
    res.innerPaths = innerPaths;
    res.runtimeSeconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() *
        1e-3;
    return res;
}

} // namespace QuantExt
