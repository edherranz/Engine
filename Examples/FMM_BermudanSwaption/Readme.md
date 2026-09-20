# FMM Bermudan swaption

A 5-year annual-exercise USD-SOFR payer swaption on the pinned `Examples/Products` market
(as of 2025-02-10), priced twice on the same trade and calibration basket (coterminal
deal-strike swaptions, bootstrap):

- `Input/ore.xml`: the generalized Forward Market Model (`Model` FMM, `Engine` LSM), three
  factors with the `RhoInf` / `Beta` correlation, displaced-diffusion volatilities with the
  `1/tau` shift, Longstaff-Schwartz lower bound plus the Andersen-Broadie dual bound. Writes
  `Output/npv.csv` and `Output/additional_results.csv` (`fmm*` keys: calibration table,
  exercise probabilities, lower / upper bound, duality gap, grid diagnostics).
- `Input/ore_lgm.xml`: ORE's LGM grid engine, output in `Output_LGM`.

Market inputs (quotes, curve configuration, conventions) are those of `Examples/Products`;
`Input/todaysmarket.xml` only selects the USD-SOFR curve, the SOFR swap index family and the
USD-SOFR normal swaption smile. See `ore-fmm/docs/A5_PLAN.md` for the configuration contract.

Run `python run.py`.
