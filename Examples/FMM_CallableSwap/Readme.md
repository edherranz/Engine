# FMM callable swap with a notice period

A 5-year USD-SOFR payer swap (pay 4.25% fixed, receive SOFR compounded in arrears) that the
holder may cancel on every anniversary from year one, deciding 30 calendar days before the
anniversary (`NoticePeriod` 30D, US calendar, Preceding). ORE represents the trade as the swap
plus a Bermudan swaption on the reverse swap whose exercise dates are the notice dates and whose
settlement dates are the anniversaries; the FMM engine decides at the notice grid date and enters
the flows after the settlement grid date (the A4 notice-period rights, adapted payoff in the dual
bound). `Input/ore.xml` prices with the FMM (LSM plus dual bound, coterminal deal-strike
bootstrap), `Input/ore_lgm.xml` with ORE's LGM grid engine. Market inputs as in
`FMM_BermudanSwaption`.

Run `python run.py`.
