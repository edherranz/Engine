# Callable bond: the LGM Grid engine drops the coupon paid on the first event date

Minimal reproduction on the stock example market (`Examples/Input`, as of 2016-02-05). The
portfolio holds a 5-year EUR 2% annual bond issued on the valuation date, callable at par on every
coupon date from year one (no credit risk), and the same bond as a plain `Bond` trade whose NPV is
the straight-bond reference.

Three runs price the callable bond with different engines and report `strippedBondNpv`
(the straight bond as seen by the callable engine) and `callPutValue`:

- `Input/ore_grid.xml`: `LGM` / `Grid` (`NumericLgmCallableBondEngine` on the event-time grid of
  the convolution solver, no intermediate time steps)
- `Input/ore_fd.xml`: `LGM` / `FD` (same engine class, finite-difference solver, 24 time steps per
  year)
- `Input/ore_mc.xml`: `CrossAssetModel` / `MC` (Monte Carlo, Longstaff-Schwartz)

Expected: `strippedBondNpv` equals the plain bond's NPV for every engine. Observed: the Grid run is
short by exactly the discounted first coupon (the coupon paid on the first call date), and its
callable NPV by the same amount; the embedded call value (`callPutValue`) is unaffected.

Mechanism (`QuantExt/qle/pricingengines/numericlgmcallablebondengine.cpp`, step 9.2): at its own
pay time a coupon has `couponRatio == 0` and is excluded from the underlying, on the assumption that
an earlier rollback step collects it. The loop runs from the last grid time down to the first event
and never processes the reference date, so on an event-only time grid there is no earlier step for
a coupon paid at the first event after the valuation date; its status stays `Open` and it is in
neither the cache nor the provisional sum when the result is assembled. Finite-difference runs are
unaffected because their intermediate time steps collect the coupon before its pay time. A candidate
fix that adds the time-zero value of such coupons before the result is assembled accompanies this
folder.

The Monte Carlo engine needs a credit curve, so `portfolio_credit.xml` and the `ore_*_credit.xml`
run files book the same bond with the stock market's `BOND_YIELD_EUR_OVER_OIS` curve for all three
engines; only the Grid run is short of the plain bond there as well (FD and MC agree with it and
with each other).

Run `python run.py` (or `ore Input/ore_grid.xml` etc. from this folder).
