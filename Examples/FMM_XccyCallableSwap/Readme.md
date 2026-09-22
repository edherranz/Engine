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

USD discounted and projected on SOFR; EUR discounted on ESTR, EURIBOR forwarding curves behind
the EUR swaption surface (ORE's EUR swap index family is EURIBOR based, the standard EUR
floating coupon, discounted on the overnight curve); EUR/USD spot and volatility surface; three
synthetic rate/FX correlations (`Input/curveconfig_xccy.xml`, `Input/marketdata_xccy.csv`:
rate/rate 0.5, FX/USD rate -0.2, FX/EUR rate 0.2, plus aliases for the label ORE uses when the
USD leg is fixed). Both curve configurations and both market data files are passed as
comma-separated lists in the run files.

## Runs

- `Input/ore_cam.xml`: stock ORE comparator, `BermudanSwaption_XCcy` on the CrossAssetModel
  (`Model` LGM, `Engine` MC in `Input/pricingengine_cam.xml`: one-factor LGM per currency
  bootstrapped to coterminal ATM swaptions, lognormal FX bootstrapped to ATM-forward EUR/USD
  options on the exercise dates, Longstaff-Schwartz with 50,000 training and 50,000 pricing
  paths, regression order 4).
- `Input/ore_cam_reduced.xml`: the telescoped representation on the same model.
- FMM runs (Hagan reduction on an EUR-ESTR FMM, joint USD/EUR FMM with FX) are added by A9.1 and
  A9.5; see `docs/A9_XCCY_TEMPLATE_PLAN.md`.

Run `python run.py` inside the environment that holds the fork Python package (the runs go
through `ORE.OREApp`, no executable).
