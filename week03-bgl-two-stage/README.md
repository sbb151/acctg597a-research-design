# Week 3: Two Stages, One Regression

**Paper:** Beaver, Griffin, and Landsman (1982), "The Incremental Information Content of Replacement Cost Earnings," *Journal of Accounting and Economics* 4(1): 15–39.

## The design flaw

BGL ask whether replacement cost earnings (PRE, from SEC ASR 190 disclosures) explain cross-sectional stock returns beyond historical cost earnings (HC). The two variables are 0.84 correlated. To handle the collinearity, BGL introduce a two-stage residual regression run in both directions: regress PRE on HC, keep the residual Z, then regress returns on HC and Z (Panel A); then reverse the roles (Panel B). The paper advertises the design as one that "permits a determination of the incremental explanatory power of collinear variables" (p. 15), and it was cited on that basis (Beaver and Landsman 1983; Lev and Ohlson 1982).

The problems, in the order the simulation demonstrates them:

1. **Correlation tables hold nothing fixed.** From BGL's own Table 3 (corr(HC, PRE) = 0.84, corr(HC, R) = 0.37, corr(PRE, R) = 0.30), the partial correlation of PRE with returns holding HC fixed is −0.02. The pairwise 0.30 answers no ceteris paribus question, and the Table 4 result is already implied by the three numbers in Table 3.
2. **Both panels of Table 4 are one regression read twice.** Substitute the first stage PRE = a + c·HC + Z into R = b0 + b1·HC + b2·PRE + e and the second stage becomes R = (b0 + b2·a) + (b1 + c·b2)·HC + b2·Z + e. The Panel (A) residual coefficient is b2, the ordinary multiple-regression coefficient on PRE, with the same t-statistic; the Panel (A) HC coefficient is the bivariate slope of R on HC. Panel (B) returns b1 and the other bivariate slope. BGL's footnote 8 concedes the equivalence. Recreated from their Table 4, the one regression they never print is R = α + 0.39·HC (t = 4.2) − 0.02·PRE (t = −0.4), R² = 0.14.
3. **Orthogonalizing does not reduce the standard error.** The two-stage residual coefficient equals the one-regression coefficient in every sample, so its sampling distribution is identical and follows σ/(√N·sd(PRE)·√(1 − ρ²)). At ρ = 0.84 the variance inflation factor is 3.4 in either procedure. Collinearity is a property of the data, not the estimator (Christie, Kennelley, King, and Schaefer 1984).
4. **Simpson's paradox.** BGL report 1977, 1978, and pooled slopes of 0.36, 0.29, and 0.32 and treat the agreement as reassurance. With the same firm-level draws in both years and only the year means shifted, identical year slopes of 0.34 pool to −0.20, and zero year slopes pool to +0.34. The pooled slope is a sums-of-squares weighted average of the year slopes plus a between-year term; the informative comparison is a test of slope equality, which BGL do not report.

## Files

- `bgl_two_stage.do` — self-contained simulation, base Stata only (no user-written packages, no data downloads). Runtime about two minutes. Run from this folder: `do bgl_two_stage.do`.
- `figures/fig1_pairwise_partial.png` — pairwise versus partial correlations, recreated from BGL Table 3.
- `figures/fig2_added_variable.png` — raw scatter of returns on PRE versus the Frisch-Waugh-Lovell added-variable plot holding HC fixed, simulated data calibrated to BGL Tables 2–4.
- `figures/fig3_collinearity_se.png` — Monte Carlo sampling standard deviation of the incremental coefficient under one-regression and two-stage estimation, against the correlation between HC and PRE.
- `figures/fig4_simpson.png` — Simpson's paradox: identical year slopes with a reversed pooled slope, and zero year slopes with a positive pooled slope.
- `bgl_two_stage.log` — full Stata log from the run that produced the committed figures, including the identity checks (all zero to six decimals).

## Calibration

Simulation parameters come from BGL Tables 2–4 for 1977: sd(HC) = 0.216, sd(PRE) = 0.377, corr(HC, PRE) = 0.84, one-regression coefficients 0.39 on HC and −0.02 on PRE, R² = 0.14, 303 firms. The seed for Module 2 was chosen so that a single draw approximately reproduces BGL's correlations, R², and t-statistics (Panel A: 0.38 (7.1) and −0.02 (−0.5) versus BGL's 0.36 (7.1) and −0.02 (−0.4); Panel B: 0.17 (5.6) and 0.42 (4.3) versus 0.16 (5.7) and 0.39 (4.2)).

## The citation trail

- BGL (1982) introduce the two-stage design and, in footnote 8, concede that it "is equivalent to an F test on adding PRE to the regression of return (R) on historical cost earnings (HC)."
- Beaver and Landsman (1983, p. 66) and Lev and Ohlson (1982, p. 265) cite the design as circumventing collinearity.
- Christie, Kennelley, King, and Schaefer (1984, *JAE*) prove the equivalence, show the raw-variable coefficient is a bivariate slope that absorbs all shared explanatory power, and argue that BGL's null on PRE may reflect collinearity rather than irrelevance.
- Jennings (1990, *TAR*) reframes "incremental information content": across the three possible two-regressor specifications under an accounting identity, only three of the six coefficients are unique, so the meaning of an incremental test depends on the conditioning variable.

## Key readings

| Role | Paper |
|---|---|
| The design | Beaver, Griffin, and Landsman (1982, *JAE*) |
| The endorsements | Lev and Ohlson (1982, *JAR* Supplement); Beaver and Landsman (1983, FASB) |
| The correction | Christie, Kennelley, King, and Schaefer (1984, *JAE*) |
| The reframing | Jennings (1990, *TAR*) |
| Interpreting the null | Burks (2023, SSRN), effect size confidence |
