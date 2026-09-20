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

/*! \file fmmirmodel.hpp
    \brief IrModel adapter over the ForwardMarketModel so that OREData's IrModelBuilder machinery
           (lazy recalibration, swaption baskets, parameter reset) applies to the FMM without
           CrossAssetModel integration (ore-fmm docs/A5_PLAN.md section 2). State-process style
           access (CAM simulation) is deliberately unsupported in this release and fails fast.
    \ingroup models
*/

#ifndef quantext_fmm_ir_model_hpp
#define quantext_fmm_ir_model_hpp

#include <qle/models/forwardmarketmodel.hpp>
#include <qle/models/irmodel.hpp>

namespace QuantExt {

class FmmIrModel : public IrModel {
public:
    explicit FmmIrModel(const QuantLib::ext::shared_ptr<ForwardMarketModel>& model);

    const QuantLib::ext::shared_ptr<ForwardMarketModel>& fmm() const { return model_; }
    const QuantLib::ext::shared_ptr<FmmParametrization>& parametrization() const { return model_->parametrization(); }

    //! \name IrModel interface
    //@{
    Measure measure() const override { return Measure::BA; }
    const QuantLib::ext::shared_ptr<Parametrization> parametrizationBase() const override {
        return model_->parametrization();
    }
    Handle<YieldTermStructure> termStructure() const override { return model_->parametrization()->termStructure(); }
    Size n() const override { return model_->parametrization()->numberOfRates(); }
    Size m() const override { return model_->parametrization()->factors(); }
    Size n_aux() const override { return 0; }
    Size m_aux() const override { return 0; }
    QuantLib::ext::shared_ptr<StochasticProcess> stateProcess() const override;
    Real discountBond(const Time t, const Time T, const Array& x,
                      const Handle<YieldTermStructure>& discountCurve = Handle<YieldTermStructure>()) const override;
    Real numeraire(const Time t, const Array& x,
                   const Handle<YieldTermStructure>& discountCurve = Handle<YieldTermStructure>()) const override;
    Real shortRate(const Time t, const Array& x,
                   const Handle<YieldTermStructure>& discountCurve = Handle<YieldTermStructure>()) const override;
    Array marginalStep(const Time t0, const Array& x0, const Time dt, const Array& dw) const override;
    QuantLib::ext::shared_ptr<IrModel> clone() const override;
    //@}

protected:
    void generateArguments() override;

private:
    QuantLib::ext::shared_ptr<ForwardMarketModel> model_;
};

} // namespace QuantExt

#endif
