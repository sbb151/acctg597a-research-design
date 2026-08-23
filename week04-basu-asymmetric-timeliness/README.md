# Week 4: The Basu (1997) Legacy — Spurious Asymmetric Timeliness

**Paper:** Basu (1997), "The Conservatism Principle and the Asymmetric Timeliness of Earnings," *Journal of Accounting and Economics* 24(1): 3–37.

Companion to the Week 1 module (`week01-reverse-regression/`), which shows that the reverse regression of earnings on returns — introduced by Beaver, Lambert, and Morse (1980) and defended by Beaver, Lambert, and Ryan (1987) — recovers the forward coefficient only when R² = 1, and that grouping by the dependent variable is mechanically biased. Basu (1997), with 4,200+ citations, builds the conditional-conservatism literature on the same earnings-on-returns direction, adding a piecewise split on the sign of returns:

g = b₀ + b₁·D + b₂·G + b₃·(D·G) + e,  D = 1(G < 0),

with b₃ > 0 read as conservatism (earnings timelier for bad news).

## The design flaw

The Basu regression inherits every pathology of the reverse regression, plus one of its own: it conditions on the sign of the regressor (returns), a variable correlated with the regression error. In a simulated world with **perfectly symmetric earnings timeliness** and right-skewed non-earnings return noise (limited liability; growth-option and discount-rate news produce occasional large positive returns unrelated to current earnings), negative returns are disproportionately news-driven while large positive returns are disproportionately noise-driven, so the reverse slope cov(g,G)/var(G) is steeper on the negative side. The Basu coefficient measures the variance composition of returns, not accounting.

**Results at N = 700 (1,000 Monte Carlo samples):** mean spurious b₃ = 0.21, mean t = 5.2, **100% rejection** of the true null. Placebo with symmetric non-earnings noise: 5.4% rejection, exactly nominal size — confirming the skewness of the return noise, not accounting behavior, drives the result.

## Files

- `basu_spurious_timeliness.do` — self-contained simulation, base Stata only. Runtime ≈ 1 minute. Run from this folder: `do basu_spurious_timeliness.do`.
- `figures/fig3_spurious_basu.png` — sampling distribution of the spurious b₃ (true b₃ = 0 in both worlds; the skewed-noise distribution never touches zero).
- `figures/fig4_basu_piecewise.png` — the piecewise fit that "finds" conservatism in symmetric data.
- `basu_spurious_timeliness.log` — full Stata log from the run that produced the committed figures.

## Key readings

- Dietrich, Muller, and Riedl (2007, *RAST*), "Asymmetric Timeliness Tests of Accounting Conservatism" — the bridge to the econometrics critique (the only accounting paper citing both BLR 1987 and Goldberger 1984); simulation evidence of spurious asymmetric timeliness.
- Givoly, Hayn, and Natarajan (2007, *TAR*); Patatoukas and Thomas (2011, 2016, *TAR*) — placebo evidence.
- Ball, Kothari, and Nikolaev (2013, *TAR* and *JAR*) — the defense.
- Collins, Hribar, and Tian (2014, *JAE*) — cash-flow asymmetry: much of the measured asymmetry cannot be accrual conservatism.
- Dietrich, Muller, and Riedl (2023, *RAST*) — validity tests: designs built on asymmetric timeliness reject a true null far too often, and proposed fixes do not restore validity.
