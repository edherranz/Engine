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

#include <qle/models/fmmirmodel.hpp>

namespace QuantExt {

namespace {
[[noreturn]] void standaloneOnly(const char* what) {
    QL_FAIL("FmmIrModel::" << what
                           << ": the FMM is a standalone pricing model in this release (no CrossAssetModel "
                              "state process; see ore-fmm A9.2). Use FmmPathGenerator / FmmLsmPricer.");
}
} // namespace

FmmIrModel::FmmIrModel(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model) : model_(model) {
    QL_REQUIRE(model_, "FmmIrModel: model is null");
    const auto& p = *model_->parametrization();
    // replication mode carries no free parameters (they belong to the wrapped LGM parametrization)
    arguments_.resize(p.numberOfParameters());
    for (Size i = 0; i < p.numberOfParameters(); ++i)
        arguments_[i] = p.parameter(i);
    registerWith(p.termStructure());
}

void FmmIrModel::generateArguments() { model_->parametrization()->update(); }

QuantLib::ext::shared_ptr<StochasticProcess> FmmIrModel::stateProcess() const { standaloneOnly("stateProcess"); }

Real FmmIrModel::discountBond(const Time, const Time, const Array&, const Handle<YieldTermStructure>&) const {
    standaloneOnly("discountBond");
}

Real FmmIrModel::numeraire(const Time, const Array&, const Handle<YieldTermStructure>&) const {
    standaloneOnly("numeraire");
}

Real FmmIrModel::shortRate(const Time, const Array&, const Handle<YieldTermStructure>&) const {
    standaloneOnly("shortRate");
}

Array FmmIrModel::marginalStep(const Time, const Array&, const Time, const Array&) const {
    standaloneOnly("marginalStep");
}

QuantLib::ext::shared_ptr<IrModel> FmmIrModel::clone() const {
    // the parametrization object is shared by design (lazy recalibration relinks the model, it
    // does not copy it); a clone therefore shares parameters, as LGM::clone() does for its
    // parametrization
    return QuantLib::ext::make_shared<FmmIrModel>(QuantLib::ext::make_shared<ForwardMarketModel>(
        model_->parametrization(), model_->measure(), model_->terminalIndex()));
}

} // namespace QuantExt
