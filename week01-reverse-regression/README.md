# Week 1: The Reverse Regression Fallacy

**Paper:** Beaver, Lambert, and Ryan (1987), "The Information Content of Security Prices: A Second Look," *Journal of Accounting and Economics* 9(2): 139–157.

## The design flaw

BLR (1987) argue that regressing earnings changes on returns (the "reverse regression," g = δG + e) mitigates the errors-in-variables problem afflicting the forward return-earnings regression (G = γg + u), because reversing the regression "places the measurement error in the disturbance term" (p. 143). They also formally derive that Beaver, Lambert, and Morse's (1980) practice of grouping data into portfolios by the *dependent* variable is asymptotically equivalent to the reverse regression.

The problems, in order of appearance in the simulation:

1. **The reverse regression recovers the forward coefficient only when R² = 1.** BLR's own Property (3) is γδ = R², so the implied forward coefficient 1/δ = γ/R². At the R² ≈ 0.13 typical of return-earnings regressions, the implied coefficient overstates γ by a factor of 1/0.13 ≈ 7.7. Goldberger (1984) had already shown, three years before BLR, that the reverse regression estimates a *different parameter* — not a bias-corrected version of the same one.
2. **The "measurement error" is not measurement error.** The garbling in BLR's framework afflicts a latent construct (permanent/economic income), not observed earnings — reported earnings are what they are, measured exactly. Meanwhile returns carry abundant *non-earnings* information, so the dependent variable in the reverse regression carries its own "error," breaking the G ⊥ e independence the whole argument requires (a point BLR themselves concede on p. 151). Two economies — one with a large transitory earnings component, one with none — produce *identical* forward, reverse, and grouped regression output: δ is not identified as a structural parameter.
3. **Grouping by the dependent variable is mechanical inflation, not error diversification.** Sorting portfolios on G sorts partly on the disturbance u itself, inflating the grouped slope toward γ/R² even when there is *no measurement error anywhere in the system*. The simulation reproduces BLR's Table 1 escalation (their 0.31 → 1.13 → 1.58 → 1.83 → 2.07 → 2.30; simulated 0.31 → 1.00 → 1.41 → 1.75 → 2.07 → 2.15) from an error-free world. Wald (1940) and Madansky (1959) require grouping to be independent of the errors; Pakes (1982) proves data-dependent grouping is asymptotically biased.
The design was institutionalized rather than corrected: Basu (1997) — 4,200+ citations — builds conditional conservatism on the same earnings-on-returns direction. The companion module `week04-basu-asymmetric-timeliness/` demonstrates the spurious asymmetric timeliness that results; it accompanies the Week 4 reading of Basu (1997).

## Files

- `reverse_regression.do` — self-contained simulation, base Stata only (no user-written packages, no data downloads). Runtime ≈ 1–2 minutes. Run from this folder: `do reverse_regression.do`.
- `figures/fig1_implied_vs_r2.png` — implied coefficient 1/δ = γ/R² explodes as R² falls; unbiased only at R² = 1.
- `figures/fig2_grouping_escalation.png` — grouping by the dependent variable reproduces BLR's Table 1 escalation with zero measurement error.
- `reverse_regression.log` — full Stata log from the run that produced the committed figures.

## Calibration

All simulation parameters come from BLR's Tables 1–2: forward γ̂ = 0.31, mean R² = 0.13, predicted grouped asymptote γ/R² = 2.38, annual cross-sections of ≈ 700 NYSE firms over 19 years (1965–1983). One note: BLR report mean δ̂ = 0.59, which is internally inconsistent with their own Property (3) (0.31 × 0.59 = 0.18 ≠ 0.13) because their reported means are averages of yearly ratios; the simulation matches γ and R², implying δ = R²/γ ≈ 0.42.

## The citation trail (the "institutionalization" evidence)

- The econometrics critique predates the accounting adoption: Conway and Roberts (1983, *JBES*) propose reverse regression as a measurement-error fix; Goldberger (1984, *JHR*; 1984, *JBES* "Redirecting Reverse Regression") and Greene (1984, *JBES*) refute it; Klepper and Leamer (1984, *Econometrica*) show even the forward/reverse "bounding" defense holds only in the classical one-regressor EIV model.
- BLR (1987) cite Conway and Roberts (1983) — but not Goldberger's refutation of it.
- Basu (1997) cites BLM (1980), not BLR (1987), and no Goldberger.
- Dietrich, Muller, and Riedl (2007, *RAST*) is the unique bridge paper citing both BLR (1987) and Goldberger (1984); the subsequent debate (Patatoukas and Thomas 2011, 2016; Ball, Kothari, and Nikolaev 2013 ×2; Collins, Hribar, and Tian 2014; Dietrich, Muller, and Riedl 2023) proceeded largely without re-engaging the econometrics literature.

## Key readings

| Role | Paper |
|---|---|
| The temptation | Conway and Roberts (1983, *JBES*) |
| The refutation | Goldberger (1984, *JHR*); Goldberger (1984, *JBES*); Greene (1984, *JBES*) |
| The bounds caveat | Klepper and Leamer (1984, *Econometrica*) |
| Grouping conditions | Wald (1940); Madansky (1959); Pakes (1982) |
| The adoption | Beaver, Lambert, and Morse (1980, *JAE*); Beaver, Lambert, and Ryan (1987, *JAE*) |
| The institutionalization | Basu (1997, *JAE*) |
| The reckoning | Dietrich, Muller, and Riedl (2007, 2023, *RAST*); Patatoukas and Thomas (2011, 2016, *TAR*) |
