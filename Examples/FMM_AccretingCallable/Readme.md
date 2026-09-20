# FMM accreting callable note

A 5-year USD note whose notional accretes 5% a year (10.0m to 12.155m); coupons of 4.5%
semiannual 30/360 on the accreted notional, callable at par of the accreted notional (dirty,
accruals included) on every coupon date from year one; the accreted amount is repaid at call or
maturity. `Input/ore.xml` prices with the FMM issuer-cancel engine, `Input/ore_lgm.xml` with
ORE's LGM FD engine. Market inputs as in `FMM_CallableBond`.

Run `python run.py`.
