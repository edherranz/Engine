# FMM callable bond

A 5-year 4.5% semiannual USD bond issued on the valuation date (2025-02-10), callable at par
(dirty, accruals included) on every coupon date from year one, valued off the USD-SOFR curve of
the pinned `Examples/Products` market without credit risk:

- `Input/ore.xml`: the generalized Forward Market Model (`Model` FMM, `Engine` LSM) in the
  issuer-cancel representation (holder value = stripped bond minus the issuer's call), with the
  underlying control variate and the Andersen-Broadie dual bound; calibrated by bootstrap to the
  coterminal ATM swaptions on the `ReferenceCalibrationGrid`. Output in `Output/`
  (`strippedBondNpv`, `callPutValue`, `fmm*` keys).
- `Input/ore_lgm.xml`: ORE's LGM finite-difference callable bond engine, output in `Output_LGM`.

The model grid is built from the reference-grid anniversaries plus the maturity, filled by the
3M `Grid` tenor; the bond's business-day-adjusted coupon and call dates map onto it within
`GridToleranceDays`. See `ore-fmm/docs/A5_PLAN.md` section 10 for the contract and the known
limitations.

Run `python run.py`.
