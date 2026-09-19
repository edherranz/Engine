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

#include <ql/math/matrixutilities/qrdecomposition.hpp>
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
            const Real flow =
                instrument_.fixedFlows[j] + instrument_.floatWeights[j] * p.tau(j) * st.R[j - 1];
            d.deflatedFlows[j] = flow / bank[j];
            total += d.deflatedFlows[j];
        }
        d.regressors.resize(nRights);
        d.deflatedFee.resize(nRights);
        d.switchedValue.resize(nRights);
        for (Size r = 0; r < nRights; ++r) {
            const auto& rt = instrument_.rights[r];
            const auto& st = path.states[rt.noticeIdx - 1]; // state at T_notice
            // realized deflated value of the flows switched by exercise (j > settleIdx)
            Real sw = 0.0;
            for (Size j = rt.settleIdx + 1; j <= last; ++j)
                sw += d.deflatedFlows[j];
            d.switchedValue[r] = sw;
            d.deflatedFee[r] =
                rt.feeFlow * model_->discountBond(st, p.rateTime(rt.settleIdx)) / bank[rt.noticeIdx];
            // regressors: par rate to maturity, front rate, deflated current mark of switched flows
            Real annuity = 0.0, mark = 0.0;
            const Real pEnd = model_->discountBond(st, p.rateTime(last));
            for (Size j = rt.noticeIdx + 1; j <= last; ++j)
                annuity += p.tau(j) * model_->discountBond(st, p.rateTime(j));
            for (Size j = rt.settleIdx + 1; j <= last; ++j) {
                const Real pj = model_->discountBond(st, p.rateTime(j));
                const Real pjm = model_->discountBond(st, p.rateTime(j - 1));
                mark += instrument_.fixedFlows[j] * pj + instrument_.floatWeights[j] * (pjm - pj);
            }
            Array x(3);
            x[0] = annuity > QL_EPSILON ? (1.0 - pEnd) / annuity : 0.0;
            x[1] = st.R[std::min<Size>(rt.noticeIdx + 1, M_) - 1];
            x[2] = mark / bank[rt.noticeIdx];
            d.regressors[r] = x;
        }
        d.deflatedFlows[0] = total; // cache the never-exercise total in slot 0
    }
}

Real FmmLsmPricer::pathValue(const PathData& d, const FmmLsmPolicy& pol, Integer* exercisedRight) const {
    const bool enter = instrument_.style == FmmCallableInstrument::Style::Enter;
    const Real total = d.deflatedFlows[0];
    if (exercisedRight)
        *exercisedRight = -1;
    for (Size r = 0; r < instrument_.rights.size(); ++r) {
        if (pol.coefficients[r].empty())
            continue;
        if (enter && d.regressors[r][2] <= 0.0)
            continue; // outside the trained domain
        const Array b = basis(d.regressors[r]);
        Real dHat = 0.0;
        for (Size k = 0; k < b.size(); ++k)
            dHat += b[k] * pol.coefficients[r][k];
        if (dHat < 0.0) {
            if (exercisedRight)
                *exercisedRight = static_cast<Integer>(r);
            return enter ? d.switchedValue[r] + d.deflatedFee[r]
                         : total - d.switchedValue[r] + d.deflatedFee[r];
        }
    }
    return enter ? 0.0 : total;
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
        // regression target: continuation minus exercise
        std::vector<Size> idx;
        for (Size n = 0; n < train.size(); ++n)
            if (!enter || train[n].regressors[rr][2] > 0.0)
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
            const Real exValue = enter ? d.switchedValue[rr] + d.deflatedFee[rr]
                                       : d.deflatedFlows[0] - d.switchedValue[rr] + d.deflatedFee[rr];
            y[q] = V[idx[q]] - exValue;
        }
        policy_.coefficients[rr] = qrSolve(A, y);
        for (const Size n : idx) {
            const PathData& d = train[n];
            const Array b = basis(d.regressors[rr]);
            Real dHat = 0.0;
            for (Size k = 0; k < nb; ++k)
                dHat += b[k] * policy_.coefficients[rr][k];
            if (dHat < 0.0)
                V[n] = enter ? d.switchedValue[rr] + d.deflatedFee[rr]
                             : d.deflatedFlows[0] - d.switchedValue[rr] + d.deflatedFee[rr];
        }
    }
    Real trainMean = 0.0;
    for (const Real v : V)
        trainMean += v;
    trainMean /= static_cast<Real>(V.size());

    // ---- independent valuation pass with the frozen policy ----
    FmmLsmResult res = valueWithPolicy(policy_, config_.valuationSeed);
    res.trainingValue = trainMean;

    // t=0 curve value of the underlying flows (Cancel style diagnostic)
    if (!enter) {
        const auto& p = *model_->parametrization();
        Real u = 0.0;
        for (Size j = 1; j <= instrument_.lastFlowIdx; ++j) {
            const Real pj = p.termStructure()->discount(p.rateTime(j));
            const Real pjm = p.termStructure()->discount(p.rateTime(j - 1));
            u += instrument_.fixedFlows[j] * pj + instrument_.floatWeights[j] * (pjm - pj);
        }
        res.underlyingValue = u;
    }
    res.runtimeSeconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() *
        1e-3;
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
    FmmLsmResult res;
    IncrementalStatistics stats;
    std::vector<IncrementalStatistics> exStats(instrument_.rights.size());
    IncrementalStatistics neverStats;
    Real timeSum = 0.0;
    Size timeCount = 0;
    const auto& p = *model_->parametrization();
    for (const auto& d : val) {
        Integer er = -1;
        stats.add(pathValue(d, pol, &er));
        for (Size r = 0; r < instrument_.rights.size(); ++r)
            exStats[r].add(er == static_cast<Integer>(r) ? 1.0 : 0.0);
        neverStats.add(er < 0 ? 1.0 : 0.0);
        if (er >= 0) {
            timeSum += p.rateTime(instrument_.rights[er].noticeIdx);
            ++timeCount;
        }
    }
    res.lowerBound = stats.mean();
    res.lowerBoundSe = stats.errorEstimate();
    for (Size r = 0; r < instrument_.rights.size(); ++r) {
        res.exerciseProbability.push_back(exStats[r].mean());
        res.exerciseProbabilitySe.push_back(exStats[r].errorEstimate());
    }
    res.noExerciseProbability = neverStats.mean();
    res.expectedExerciseTime = timeCount > 0 ? timeSum / static_cast<Real>(timeCount) : 0.0;
    return res;
}

std::pair<Real, Real> FmmLsmPricer::pairedPolicyDifference(const FmmLsmPolicy& importedPolicy,
                                                           const BigNatural seed) const {
    QL_REQUIRE(importedPolicy.coefficients.size() == instrument_.rights.size(),
               "FmmLsmPricer::pairedPolicyDifference: policy rights do not align");
    std::vector<PathData> val;
    simulate(config_.valuationPaths, seed, config_.valuationSequence, val);
    IncrementalStatistics diff;
    for (const auto& d : val)
        diff.add(pathValue(d, policy_, nullptr) - pathValue(d, importedPolicy, nullptr));
    return {diff.mean(), diff.errorEstimate()};
}

} // namespace QuantExt
