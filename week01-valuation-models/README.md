# Week 1 companion: The Valuation Model the ERC Literature Assumes

**Paper:** Easton and Zmijewski (1989), "Cross-Sectional Variation in the Stock Market Response to Accounting Earnings Announcements," *Journal of Accounting and Economics* 11(2-3): 117-141.

Companion to the Week 1 discussion of valuation models. EZ open with the general Garman-Ohlson (1980) model of valuation under uncertainty (price as a function of *all* dividend-predicting information), then pivot to an earnings-capitalization implementation whose assumptions are delegated entirely to a footnote citation of Ohlson (1983). Earnings capitalization (MV = E/r) and book-value capitalization (MV = BV) are valuation models that hold in a world of certainty (or knife-edge cases: earnings a sufficient statistic for payoffs; accounting that marks book to value).

## The figure

`valuation_erc_figure.do` (base Stata, seconds to run) plots the PVED-consistent theoretical ERC as a function of earnings persistence ω at discount rate r:

ERC(ω, r) = 1 + ω/(1 + r − ω),

which nests the polar cases: purely transitory earnings (ω = 0) give ERC = 1, and permanent earnings (ω = 1) give the capitalization multiple ERC = (1+r)/r ≈ 11 at r = 10%. Annual earnings are close to a random walk, so taking the capitalization model seriously predicts ERCs near 11; Easton and Zmijewski's own estimates are 1.65 (two-day window) and 2.53 (forecast-date window), and typical estimated ERCs sit between 1 and 3.

## Key readings

- Garman and Ohlson (1980, *JAR*) — the general arbitrage-free model (84 citations).
- Ohlson (1983, *JAR*) — conditions for price to equal capitalized earnings under uncertainty; the paper EZ's footnote 8 leans on (53 citations).
- Miller and Modigliani (1961, *J. Business*) — the capitalization representation and its certainty/no-growth conditions.
- Ohlson (1991, 1995, *CAR*) — the bridge: clean surplus + linear information dynamics make price a weighted average of capitalized earnings and book value, with MV = E/r and MV = BV as the ω → 1 and ω → 0 endpoints.
- Kothari and Zimmerman (1995, *JAE*) — price and return models: what each specification assumes.
- Holthausen and Watts (2001) vs. Barth, Beaver, and Landsman (2001, both *JAE*) — the debate over running regressions whose valuation theory is assumed rather than developed.
