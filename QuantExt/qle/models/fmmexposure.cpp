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

#include <qle/models/fmmexposure.hpp>

#include <ql/math/matrixutilities/qrdecomposition.hpp>
#include <ql/math/matrixutilities/choleskydecomposition.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>

namespace QuantExt {

namespace {
const Time tEps = 1e-10;

struct MeanAcc {
    Real sum = 0.0, sumSq = 0.0;
    Size n = 0;
    void add(const Real v) {
        sum += v;
        sumSq += v * v;
        ++n;
    }
    Real mean() const { return n > 0 ? sum / static_cast<Real>(n) : Null<Real>(); }
    Real se() const {
        if (n < 2)
            return Null<Real>();
        const Real m = mean();
        const Real var = std::max(0.0, (sumSq / static_cast<Real>(n) - m * m) * static_cast<Real>(n) /
                                           static_cast<Real>(n - 1));
        return std::sqrt(var / static_cast<Real>(n));
    }
};
} // namespace

struct FmmExposureEngine::TradeData {
    Size idx = 0;
    const FmmExposureTrade* trade = nullptr;
    const FmmLsmPricer* helper = nullptr;
};

FmmExposureEngine::FmmExposureEngine(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model,
                                     const std::vector<FmmExposureTrade>& trades, const std::vector<Time>& exposureTimes,
                                     const FmmExposureConfig& config)
    : model_(model), trades_(trades), exposureTimes_(exposureTimes), config_(config) {
    QL_REQUIRE(model_, "FmmExposureEngine: model is null");
    QL_REQUIRE(!trades_.empty(), "FmmExposureEngine: no trades");
    QL_REQUIRE(!exposureTimes_.empty(), "FmmExposureEngine: no exposure times");
    QL_REQUIRE(config_.trainingPaths > 0 && config_.exposurePaths > 0, "FmmExposureEngine: path counts must be positive");
    QL_REQUIRE(config_.trainingSeed != config_.exposureSeed || config_.trainingSequence != config_.exposureSequence,
               "FmmExposureEngine: training and exposure paths must be independent (different seeds)");
    const auto& p = *model_->parametrization();
    const Size M = p.numberOfRates();
    for (Size i = 0; i < exposureTimes_.size(); ++i) {
        QL_REQUIRE(exposureTimes_[i] > (i == 0 ? 0.0 : exposureTimes_[i - 1]) + tEps,
                   "FmmExposureEngine: exposure times must be positive and strictly increasing");
        QL_REQUIRE(exposureTimes_[i] <= p.rateTime(M) + tEps,
                   "FmmExposureEngine: exposure time " << exposureTimes_[i] << " beyond the model grid ("
                                                        << p.rateTime(M) << ")");
    }
    FmmLsmConfig helperConfig;
    helperConfig.basisOrder = config_.basisOrder;
    for (const auto& t : trades_) {
        t.instrument.validate(M, false); // vanilla structures (no rights) are admitted
        QL_REQUIRE(t.instrument.lastFlowIdx >= 1, "FmmExposureEngine: trade " << t.id << " has no flows");
        if (!t.instrument.rights.empty())
            QL_REQUIRE(t.policy.coefficients.size() == t.instrument.rights.size(),
                       "FmmExposureEngine: trade " << t.id << ": policy has " << t.policy.coefficients.size()
                                                   << " rights, instrument " << t.instrument.rights.size());
        lastFlow_ = std::max(lastFlow_, t.instrument.lastFlowIdx);
        helpers_.push_back(QuantLib::ext::make_shared<FmmLsmPricer>(model_, t.instrument, helperConfig));
    }
    // simulation time line: the grid times up to the last flow (realised flows, decisions) and the
    // exposure times (states for the values), merged
    std::vector<std::pair<Time, int>> all; // (time, kind): kind 0 grid index j encoded below, 1 exposure
    for (Size j = 1; j <= lastFlow_; ++j)
        all.emplace_back(p.rateTime(j), static_cast<int>(j));
    for (Size d = 0; d < exposureTimes_.size(); ++d)
        all.emplace_back(exposureTimes_[d], -static_cast<int>(d) - 1);
    std::sort(all.begin(), all.end(), [](const std::pair<Time, int>& x, const std::pair<Time, int>& y) {
        return x.first < y.first;
    });
    gridPos_.assign(lastFlow_ + 1, Null<Size>());
    expPos_.assign(exposureTimes_.size(), Null<Size>());
    for (const auto& e : all) {
        if (userTimes_.empty() || e.first > userTimes_.back() + tEps)
            userTimes_.push_back(e.first);
        const Size pos = userTimes_.size() - 1;
        if (e.second > 0)
            gridPos_[static_cast<Size>(e.second)] = pos;
        else
            expPos_[static_cast<Size>(-e.second - 1)] = pos;
    }
}

bool FmmExposureEngine::isDue(const Time T, const Time t) const {
    return T > t + tEps || (config_.includeFlowsOnExposureDate && std::fabs(T - t) <= tEps);
}

Real FmmExposureEngine::markValue(const ForwardMarketModel& model, const FmmCallableInstrument& inst,
                                  const ForwardMarketModel::State& st, const Size from, const Size to,
                                  const bool includeFlowsAtT) {
    const auto& p = *model.parametrization();
    const Time t = st.t;
    const Size a = st.period; // accruing period: rates k < a are fixed, R_a(t) accrues
    QL_REQUIRE(a >= 1 && a <= p.numberOfRates(), "FmmExposureEngine::markValue: state beyond the grid");
    auto due = [&](const Size j) {
        const Time T = p.rateTime(j);
        return T > t + tEps || (includeFlowsAtT && std::fabs(T - t) <= tEps);
    };
    // (extended) bond P(t, T_j): the model bond for T_j >= t, and B(t) / B(T_{a-1}) for the start of
    // the accruing period (the extended zero-coupon bond of the generalized FMM, through which the
    // accruing rate R_a(t) reads 1 + tau_a R_a(t) = P(t, T_{a-1}) / P(t, T_a))
    auto bond = [&](const Size j) -> Real {
        if (j >= a)
            return model.discountBond(st, p.rateTime(j));
        QL_REQUIRE(j + 1 == a, "FmmExposureEngine::markValue: bond to a fixed period start " << j << " requested");
        return model.discountBond(st, p.rateTime(a)) * (1.0 + p.tau(a) * st.R[a - 1]);
    };
    auto spread = [&](const Size j) { return std::exp(-inst.issuerSpread * (p.rateTime(j) - t)); };
    Real v = 0.0;
    for (Size j = std::max(from + 1, a); j <= to; ++j) {
        if (!due(j))
            continue;
        const Real pj = bond(j), pjm = bond(j - 1);
        v += (inst.fixedFlows[j] * pj + inst.floatWeights[j] * (pjm - pj)) * spread(j);
    }
    for (const auto& f : inst.compoundedFloats) {
        if (f.payIdx <= from || f.payIdx > to || !due(f.payIdx))
            continue;
        // fully fixed periods k < a, the accruing period and the rest through the extended bond
        Real fixedFactor = 1.0;
        Size m = f.startIdx;
        for (Size k = f.startIdx + 1; k <= std::min(f.endIdx, a - 1); ++k) {
            fixedFactor *= 1.0 + p.tau(k) * st.R[k - 1];
            m = k;
        }
        const Real fwdFactor = m < f.endIdx ? bond(m) / bond(f.endIdx) : 1.0;
        v += (f.weight * (fixedFactor * fwdFactor - 1.0) + f.spreadAmount) * bond(f.payIdx) * spread(f.payIdx);
    }
    return v;
}

Real FmmExposureEngine::feeValue(const ForwardMarketModel& model, const FmmCallableInstrument& inst,
                                 const ForwardMarketModel::State& st, const Size r, const bool includeFlowsAtT) {
    const auto& p = *model.parametrization();
    const auto& rt = inst.rights[r];
    const Time T = p.rateTime(rt.settleIdx), t = st.t;
    if (!(T > t + tEps || (includeFlowsAtT && std::fabs(T - t) <= tEps)))
        return 0.0;
    return rt.feeFlow * model.discountBond(st, T) * std::exp(-inst.issuerSpread * (T - t));
}

Array FmmExposureEngine::exposureBasis(const Array& x) const {
    if (config_.basisOrder == 1)
        return Array{1.0, x[1], x[2]};
    return Array{1.0, x[1], x[2], x[1] * x[1], x[2] * x[2], x[1] * x[2]};
}

Real FmmExposureEngine::Fit::evaluate(const Array& b) const {
    Real v = yMean;
    if (!active)
        return v;
    for (Size k = 1; k < b.size(); ++k)
        v += beta[k - 1] * (b[k] - mean[k - 1]) / scale[k - 1];
    return v;
}

Real FmmExposureEngine::exactPart(const FmmCallableInstrument& inst, const Array& x) {
    // Cancel style: the remaining underlying flows are held whatever the future decisions, their
    // value is the third regressor (exact on the state's curve) and only the cancellation option is
    // regressed; Enter style: nothing is held before exercise
    return inst.style == FmmCallableInstrument::Style::Cancel ? x[2] : 0.0;
}

bool FmmExposureEngine::decisionPending(const FmmCallableInstrument& inst, const Integer exercised, const Time t) const {
    if (inst.rights.empty())
        return false;
    const auto& p = *model_->parametrization();
    const Size r = exercised >= 0 ? static_cast<Size>(exercised) : inst.rights.size() - 1;
    return p.rateTime(inst.rights[r].noticeIdx) > t + tEps;
}

Real FmmExposureEngine::committedValue(const TradeData& td, const PathTrade& pt,
                                       const ForwardMarketModel::State& st) const {
    const auto& inst = td.trade->instrument;
    const bool enter = inst.style == FmmCallableInstrument::Style::Enter;
    const bool incl = config_.includeFlowsOnExposureDate;
    if (pt.exercised < 0)
        return enter ? 0.0 : markValue(*model_, inst, st, 0, inst.lastFlowIdx, incl);
    const Size r = static_cast<Size>(pt.exercised);
    const Size settle = inst.rights[r].settleIdx;
    const Real fee = feeValue(*model_, inst, st, r, incl);
    return enter ? markValue(*model_, inst, st, settle, inst.lastFlowIdx, incl) + fee
                 : markValue(*model_, inst, st, 0, settle, incl) + fee;
}

Array FmmExposureEngine::regressors(const TradeData& td, const ForwardMarketModel::State& st) const {
    const auto& p = *model_->parametrization();
    const auto& inst = td.trade->instrument;
    const Size a = st.period, last = inst.lastFlowIdx;
    Array x(3, 0.0);
    if (a > last)
        return x;
    Real annuity = 0.0;
    for (Size j = a + 1; j <= last; ++j)
        annuity += p.tau(j) * model_->discountBond(st, p.rateTime(j));
    const Real pa = model_->discountBond(st, p.rateTime(a)), pEnd = model_->discountBond(st, p.rateTime(last));
    x[0] = annuity > QL_EPSILON ? (pa - pEnd) / annuity : 0.0;
    x[1] = st.R[a - 1];
    x[2] = markValue(*model_, inst, st, 0, last, config_.includeFlowsOnExposureDate);
    return x;
}

void FmmExposureEngine::processPath(const FmmPathGenerator::Path& path, std::vector<PathTrade>& out) const {
    const auto& p = *model_->parametrization();
    std::vector<Real> bank(lastFlow_ + 1, 1.0);
    for (Size j = 1; j <= lastFlow_; ++j)
        bank[j] = model_->bankAccount(path.states[gridPos_[j]]);
    out.assign(trades_.size(), PathTrade());
    for (Size tr = 0; tr < trades_.size(); ++tr) {
        const auto& inst = trades_[tr].instrument;
        const auto& helper = *helpers_[tr];
        PathTrade& pt = out[tr];
        const Size last = inst.lastFlowIdx;
        std::vector<Real> flows(last + 1, 0.0);
        for (Size j = 1; j <= last; ++j)
            flows[j] = helper.realizedFlow(j, path.states[gridPos_[j]]) / bank[j] *
                       std::exp(-inst.issuerSpread * p.rateTime(j));
        // the policy's decision on this path: first right at which it exercises
        pt.exercised = -1;
        for (Size r = 0; r < inst.rights.size(); ++r) {
            const Size n = inst.rights[r].noticeIdx;
            const Array x = helper.regressorsAt(path.states[gridPos_[n]], r, bank[n]);
            if (helper.exerciseDecision(x, r, trades_[tr].policy)) {
                pt.exercised = static_cast<Integer>(r);
                break;
            }
        }
        const bool enter = inst.style == FmmCallableInstrument::Style::Enter;
        pt.policyFlows.assign(last + 1, 0.0);
        if (pt.exercised >= 0) {
            const auto& rt = inst.rights[static_cast<Size>(pt.exercised)];
            for (Size j = 1; j <= last; ++j)
                if (enter ? j > rt.settleIdx : j <= rt.settleIdx)
                    pt.policyFlows[j] = flows[j];
            pt.policyFlows[rt.settleIdx] +=
                rt.feeFlow / bank[rt.settleIdx] * std::exp(-inst.issuerSpread * p.rateTime(rt.settleIdx));
        } else if (!enter) {
            pt.policyFlows = flows;
        }
        pt.total = 0.0;
        for (Size j = 1; j <= last; ++j)
            pt.total += pt.policyFlows[j];
    }
}

FmmExposureResult FmmExposureEngine::calculate() const {
    const auto start = std::chrono::steady_clock::now();
    const auto& p = *model_->parametrization();
    const Size ntr = trades_.size(), nd = exposureTimes_.size();
    std::vector<TradeData> tds(ntr);
    for (Size tr = 0; tr < ntr; ++tr) {
        tds[tr].idx = tr;
        tds[tr].trade = &trades_[tr];
        tds[tr].helper = helpers_[tr].get();
    }
    const Size nb = exposureBasis(Array(3, 0.0)).size();

    // ---- training pass: conditional-value regressions per trade and exposure date ----
    std::vector<std::vector<std::vector<Real>>> trainA(ntr, std::vector<std::vector<Real>>(nd));
    std::vector<std::vector<std::vector<Real>>> trainY(ntr, std::vector<std::vector<Real>>(nd));
    {
        FmmPathGenerator gen(model_, userTimes_, config_.trainingSequence, config_.trainingSeed);
        std::vector<PathTrade> pts;
        for (Size n = 0; n < config_.trainingPaths; ++n) {
            const auto path = gen.next();
            processPath(path, pts);
            for (Size d = 0; d < nd; ++d) {
                const auto& st = path.states[expPos_[d]];
                const Time t = st.t;
                Real bank = Null<Real>();
                for (Size tr = 0; tr < ntr; ++tr) {
                    const auto& inst = trades_[tr].instrument;
                    if (!decisionPending(inst, pts[tr].exercised, t))
                        continue;
                    if (bank == Null<Real>())
                        bank = model_->bankAccount(st);
                    Real future = 0.0;
                    for (Size j = 1; j <= inst.lastFlowIdx; ++j)
                        if (isDue(p.rateTime(j), t))
                            future += pts[tr].policyFlows[j];
                    const Real y = future * bank * std::exp(inst.issuerSpread * t);
                    const Array x = regressors(tds[tr], st);
                    const Array b = exposureBasis(x);
                    trainA[tr][d].insert(trainA[tr][d].end(), b.begin(), b.end());
                    trainY[tr][d].push_back(y - exactPart(inst, x));
                }
            }
        }
    }
    FmmExposureResult res;
    res.times = exposureTimes_;
    res.paths = config_.exposurePaths;
    res.diagnostics.assign(ntr, std::vector<FmmExposureDiagnostics>(nd));
    std::vector<std::vector<Fit>> fits(ntr, std::vector<Fit>(nd));
    for (Size tr = 0; tr < ntr; ++tr) {
        for (Size d = 0; d < nd; ++d) {
            auto& diag = res.diagnostics[tr][d];
            auto& fit = fits[tr][d];
            const Size rows = trainY[tr][d].size();
            diag.regressionPaths = rows;
            if (rows == 0)
                continue;
            Real yMean = 0.0;
            for (const Real y : trainY[tr][d])
                yMean += y;
            yMean /= static_cast<Real>(rows);
            fit.yMean = yMean;
            if (rows < std::max<Size>(config_.minRegressionPaths, 2 * nb))
                continue; // too few alive paths: the mean of their values
            // centred and scaled columns (the constant is the intercept yMean), normal equations
            // with a tiny ridge for numerical safety only
            const Size m = nb - 1;
            fit.mean = Array(m, 0.0);
            fit.scale = Array(m, 0.0);
            const auto& A = trainA[tr][d];
            for (Size q = 0; q < rows; ++q)
                for (Size k = 0; k < m; ++k)
                    fit.mean[k] += A[q * nb + k + 1];
            fit.mean /= static_cast<Real>(rows);
            for (Size q = 0; q < rows; ++q)
                for (Size k = 0; k < m; ++k) {
                    const Real c = A[q * nb + k + 1] - fit.mean[k];
                    fit.scale[k] += c * c;
                }
            for (Size k = 0; k < m; ++k)
                fit.scale[k] = std::sqrt(fit.scale[k] / static_cast<Real>(rows));
            for (Size k = 0; k < m; ++k)
                if (fit.scale[k] < 1e-14)
                    fit.scale[k] = 1.0; // a degenerate column contributes nothing
            Matrix G(m, m, 0.0);
            Array g(m, 0.0);
            for (Size q = 0; q < rows; ++q) {
                Array z(m);
                for (Size k = 0; k < m; ++k)
                    z[k] = (A[q * nb + k + 1] - fit.mean[k]) / fit.scale[k];
                const Real yc = trainY[tr][d][q] - yMean;
                for (Size k = 0; k < m; ++k) {
                    g[k] += z[k] * yc;
                    for (Size l = 0; l < m; ++l)
                        G[k][l] += z[k] * z[l];
                }
            }
            for (Size k = 0; k < m; ++k)
                G[k][k] += 1e-8 * static_cast<Real>(rows);
            fit.beta = CholeskySolveFor(CholeskyDecomposition(G), g);
            fit.active = true;
            Real ssRes = 0.0, ssTot = 0.0;
            for (Size q = 0; q < rows; ++q) {
                Array b(nb);
                for (Size k = 0; k < nb; ++k)
                    b[k] = A[q * nb + k];
                const Real e = trainY[tr][d][q] - fit.evaluate(b);
                ssRes += e * e;
                ssTot += (trainY[tr][d][q] - yMean) * (trainY[tr][d][q] - yMean);
            }
            diag.regressionR2 = ssTot > 0.0 ? 1.0 - ssRes / ssTot : Null<Real>();
            diag.regressionResidualSd = std::sqrt(ssRes / static_cast<Real>(rows));
        }
    }
    trainA.clear();
    trainY.clear();

    // ---- exposure pass ----
    const Size paths = config_.exposurePaths;
    res.numeraire.assign(nd * paths, 0.0);
    res.deflatedValues.assign(ntr, std::vector<Real>(nd * paths, 0.0));
    std::vector<MeanAcc> t0(ntr);
    std::vector<std::vector<MeanAcc>> bias(ntr, std::vector<MeanAcc>(nd)), mart(ntr, std::vector<MeanAcc>(nd));
    std::vector<std::vector<Size>> alive(ntr, std::vector<Size>(nd, 0)), noticed = alive, settled = alive,
                                   exercised = alive;
    {
        FmmPathGenerator gen(model_, userTimes_, config_.exposureSequence, config_.exposureSeed);
        std::vector<PathTrade> pts;
        for (Size n = 0; n < paths; ++n) {
            const auto path = gen.next();
            processPath(path, pts);
            for (Size tr = 0; tr < ntr; ++tr)
                t0[tr].add(pts[tr].total);
            for (Size d = 0; d < nd; ++d) {
                const auto& st = path.states[expPos_[d]];
                const Time t = st.t;
                const Real bank = model_->bankAccount(st);
                res.numeraire[d * paths + n] = bank;
                for (Size tr = 0; tr < ntr; ++tr) {
                    const auto& inst = trades_[tr].instrument;
                    const PathTrade& pt = pts[tr];
                    Real future = 0.0;
                    for (Size j = 1; j <= inst.lastFlowIdx; ++j)
                        if (isDue(p.rateTime(j), t))
                            future += pt.policyFlows[j];
                    const Real past = pt.total - future;
                    Real value;
                    if (decisionPending(inst, pt.exercised, t)) {
                        ++alive[tr][d];
                        const Array x = regressors(tds[tr], st);
                        value = exactPart(inst, x) + fits[tr][d].evaluate(exposureBasis(x));
                        bias[tr][d].add(value - future * bank * std::exp(inst.issuerSpread * t));
                    } else {
                        value = committedValue(tds[tr], pt, st);
                        if (pt.exercised >= 0) {
                            const auto& rt = inst.rights[static_cast<Size>(pt.exercised)];
                            ++exercised[tr][d];
                            if (p.rateTime(rt.settleIdx) > t + tEps)
                                ++noticed[tr][d];
                            else
                                ++settled[tr][d];
                        }
                    }
                    res.deflatedValues[tr][d * paths + n] = value / bank;
                    mart[tr][d].add(value / bank * std::exp(-inst.issuerSpread * t) + past);
                }
            }
        }
    }
    res.t0Values.resize(ntr);
    res.t0ValuesSe.resize(ntr);
    res.t0CurveValues.resize(ntr);
    for (Size tr = 0; tr < ntr; ++tr) {
        res.t0Values[tr] = t0[tr].mean();
        res.t0ValuesSe[tr] = t0[tr].se();
        res.t0CurveValues[tr] = trades_[tr].instrument.style == FmmCallableInstrument::Style::Enter
                                    ? 0.0
                                    : fmmUnderlyingCurveValue(trades_[tr].instrument, p);
        for (Size d = 0; d < nd; ++d) {
            auto& diag = res.diagnostics[tr][d];
            diag.aliveFraction = static_cast<Real>(alive[tr][d]) / static_cast<Real>(paths);
            diag.noticedFraction = static_cast<Real>(noticed[tr][d]) / static_cast<Real>(paths);
            diag.settledFraction = static_cast<Real>(settled[tr][d]) / static_cast<Real>(paths);
            diag.exercisedFraction = static_cast<Real>(exercised[tr][d]) / static_cast<Real>(paths);
            diag.outOfSampleBias = bias[tr][d].mean();
            diag.outOfSampleBiasSe = bias[tr][d].se();
            diag.martingaleCheck = mart[tr][d].mean();
            diag.martingaleCheckSe = mart[tr][d].se();
        }
    }
    res.runtimeSeconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() *
        1e-3;
    return res;
}

} // namespace QuantExt
