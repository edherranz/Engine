# FMM step-up callable note

A 5-year USD note with semiannual 30/360 coupons stepping from 4.00% to 6.00% in yearly steps,
callable at par (dirty, accruals included) on every coupon date from year one, valued off the
USD-SOFR curve without credit risk. The step-up makes the call increasingly valuable to the
issuer, so the exercise policy matters more than for the flat-coupon bond. `Input/ore.xml` prices
with the FMM issuer-cancel engine (LSM, control variate, dual bound; coterminal ATM bootstrap on
the reference calibration grid), `Input/ore_lgm.xml` with ORE's LGM FD engine. Market inputs as in
`FMM_CallableBond`.

Run `python run.py`.
