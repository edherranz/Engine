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

/*! \file fmmanalyticswaptionengine.hpp
    \brief QuantLib Swaption engine on the FMM frozen-gradient approximations (strike-flat normal
           or effective shift), used for calibration helpers and European pricing on the model grid
           (ore-fmm docs/A5_PLAN.md section 5). Physically settled, grid-aligned schedules only.
    \ingroup engines
*/

#ifndef quantext_fmm_analytic_swaption_engine_hpp
#define quantext_fmm_analytic_swaption_engine_hpp

#include <qle/models/fmmgrid.hpp>
#include <qle/models/fmmparametrization.hpp>

#include <ql/instruments/swaption.hpp>

namespace QuantExt {

class FmmAnalyticSwaptionEngine : public Swaption::engine {
public:
    FmmAnalyticSwaptionEngine(const QuantLib::ext::shared_ptr<FmmParametrization>& parametrization,
                              const QuantLib::ext::shared_ptr<FmmGrid>& grid,
                              const FmmSwaptionApproxMethod method = FmmSwaptionApproxMethod::EffectiveShift,
                              const Natural gridToleranceDays = 3);
    void calculate() const override;

    const QuantLib::ext::shared_ptr<FmmGrid>& grid() const { return grid_; }
    FmmSwaptionApproxMethod method() const { return method_; }

private:
    QuantLib::ext::shared_ptr<FmmParametrization> parametrization_;
    QuantLib::ext::shared_ptr<FmmGrid> grid_;
    FmmSwaptionApproxMethod method_;
    Natural gridToleranceDays_;
};

} // namespace QuantExt

#endif
