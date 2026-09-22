# Cancellable cross-currency swap (USD-SOFR floating against EUR fixed)

Template deal of the cross-currency milestones (A9.1 Hagan reduction, A9.5 joint two-currency
FMM): we pay USD-SOFR compounded in arrears plus 20 bp quarterly on USD 10m and receive 2.25%
annual (30/360) on EUR 9.7m, notionals exchanged at start and maturity, 5 years from 2025-02-12.
We may cancel on every anniversary from year one; on cancellation the notionals are exchanged
again and all later payments stop (Hagan, *Callable X-CCY Swap*, `docs/papers`).

Stock ORE has no cancellable cross-currency trade type. The deal is written with two stock trades
(`Input/portfolio.xml`), Hagan's decomposition `V = V_swap + V_Bermudan`:

- `XCCY_SWAP_...`: the cross-currency swap itself (`Swap` trade with legs in two currencies,
  priced by the discounting cross-currency engine);
- `XCCY_CANCEL_RIGHT_...`: a physically settled Bermudan `Swaption` into the reverse swap
  (coupons and final exchanges, no initial exchange) whose exercise fees carry the notional
  exchange on the cancellation date (`ExerciseFee currency="USD"` +10m paid, `currency="EUR"`
  -9.7m received). ORE's exercise-into logic keeps the flows after the exercise date only, so
  without the fees the option would enter the reverse coupons but never return the notionals.

`Input/portfolio_reduced.xml` is the same right with the USD floating leg telescoped away
(Hagan eqs. 1.4-1.5): under single-curve SOFR discounting the exercised USD flows (SOFR + 20 bp,
USD 10m at maturity, the USD 10m fee on the cancellation date) equal the remaining 20 bp spread
stream, so the right is a Bermudan into "receive 20 bp on USD 10m, pay 2.25% on EUR 9.7m with
the final exchange" with the EUR 9.7m rebate. Both representations must agree (consistency check
of the fee mechanism and of the telescoping identity inside ORE's engine).

## Market (`Input/todaysmarket.xml`, Products 2025-02-10 data)

USD discounted and projected on SOFR. EUR: the default pricing configuration discounts the EUR
leg on `EUR-IN-USD` (USD collateral, the market's cross-currency basis), the `inccy`
configuration on ESTR (in-currency view; the LGM / FMM calibration configuration, and the
alternative pricing view through `Markets/pricing` in the run files). EURIBOR forwarding curves
sit behind the EUR swaption surface (ORE's EUR swap index family is EURIBOR based, the standard
EUR floating coupon, discounted on the overnight curve); EUR/USD spot and volatility surface.
Correlations (`Input/curveconfig_xccy.xml`, quotes `Input/marketdata_xccy.csv`): Hagan's sizing
values, rate/rate 0.25 and FX/EUR-rate +0.25 in Hagan's convention (F^{A/B} = EUR per USD), which
is -0.25 for ORE's `FX-GENERIC-EUR-USD` (USD per EUR); FX/USD-rate 0. The low scenario
(FX/EUR-rate -0.25 in Hagan's convention) is `Input/marketdata_xccy_rhofxb_low.csv`; swap the
file name in the run file's `marketDataFile` list to use it, or edit the values. Aliases exist for
the label ORE uses when the USD leg is fixed (`USD-FedFunds`). Both curve configurations and both
market data files are passed as comma-separated lists.

## Runs

- `Input/ore_cam.xml`: stock ORE comparator, `BermudanSwaption_XCcy` on the CrossAssetModel
  (`Model` LGM, `Engine` MC in `Input/pricingengine_cam.xml`: one-factor LGM per currency
  bootstrapped to coterminal ATM swaptions on the in-currency view, lognormal FX bootstrapped to
  ATM-forward EUR/USD options on the exercise dates, Longstaff-Schwartz with 50,000 training and
  50,000 pricing paths, regression order 4).
- `Input/ore_cam_reduced.xml`: the telescoped representation on the same model; the variants
  `ore_cam_reduced_v200k4.xml` / `_v400k6.xml` with more paths and a higher regression order.
  The comparator of record is the telescoped representation at 200,000 paths averaged over ten
  seed pairs (`docs/A9_XCCY_TEMPLATE_PLAN.md` section 4).
- `Input/ore_fmm_hagan.xml`: the FMM challenger, Hagan's reduction with the improved method
  (`Input/pricingengine_fmm_hagan.xml`: `BermudanSwaption_XCcy` with `Model` FMM, `Engine` LSM,
  `XccyMethod` Improved); `Input/ore_fmm_hagan_standard.xml` the standard method. Additional
  results `fmmXccy*` report the equivalent spreads, strikes, market and effective volatilities,
  the correlations in Hagan's convention, the annuity / convexity ratios and the omitted-effect
  scale, next to the `fmm*` LSM statistics (in EUR, the modelled currency) and the calibration
  record.

Results (2026-09-22, `validation/a9_xccy_compare.py`; USD): cancellation right on the
CrossAssetModel 156,976 +/- 552 (seed mean), FMM improved method 152,864 +/- 1,349 (LSM lower
bound) to 153,929 +/- 1,366 (dual upper bound), standard method 151,751 to 152,796; the swap
+17,613. Details and the reduction diagnostics: `docs/A9_XCCY_TEMPLATE_PLAN.md` section 9.

Run `python run.py` inside the environment that holds the fork Python package (the runs go
through `ORE.OREApp`, no executable); `python run.py fmm` runs the FMM variants only.
