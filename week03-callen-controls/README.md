# Week 3: Centering, Mediators, and Residuals

**Paper:** Callen, Livnat, and Segal (2009), "The Impact of Earnings on the Pricing of Credit Default Swaps," *The Accounting Review* 84(5): 1363–1394.

## The design problems

Callen, Livnat, and Segal ask whether quarterly earnings are priced in single-name CDS premia, using levels, changes, and short-window designs on one set of regressors (structural-model determinants, the S&P rating, and restructuring-clause indicators). Three design choices in the paper are worth learning to recognize, and two of them can be demonstrated by pure simulation.

1. **Demeaning does not reduce interaction-term multicollinearity.** The paper states (p. 1373) that continuous regressors are demeaned "to reduce potential multicollinearity induced by interaction terms," citing Aiken and West (1991). Substituting x = x_c + x̄ into y = β0 + β1 x + β2 D + β3 xD + e gives y = (β0 + β1 x̄) + β1 x_c + (β2 + β3 x̄) D + β3 x_c D + e. The regressor sets span the same column space, so fitted values, R², β1, β3, and their standard errors are identical; only the intercept and the coefficient on D are relabeled. The variance inflation factor of the interaction regressor falls because the regressor's own variance falls by the same factor. Echambadi and Hess (2007, *Marketing Science*) prove the general case.
2. **The credit rating is a mediator.** S&P ratings are a function of the profitability the paper studies, so controlling for RATE (Table 3), for within-quarter rating changes (Table 4), or partitioning on the rating (Table 6) converts the hypothesized total effect of earnings on premia into a direct effect that excludes the rating channel. That direct effect is identified only if nothing unobserved moves both the rating and the premium; otherwise conditioning on the mediator induces a correlation between earnings and the shared shock (the "bad control" of Angrist and Pischke 2009; Cinelli, Forney, and Pearl 2022).
3. **Table 5 is the Beaver-Griffin-Landsman two-stage residual regression run in one direction.** The paper regresses the three-day equity return on the earnings surprise, keeps the residual RES_RET, and regresses the CDS change on the surprise and RES_RET. By the same substitution as in the Week 3 BGL module, the RES_RET coefficient equals the EQ_RET coefficient of the joint regression (−0.319 versus −0.313 in Table 5) and the EAR_SUR coefficient reverts to its bivariate slope (−1.392 versus −1.400 in column 1); adjusted R² is 0.043 in both columns. The "all-inclusive" reading assigns all shared variance to earnings by construction. This module does not simulate that identity; the `week03-bgl-two-stage` module already does.

## Files

- `callen_centering.do` — pure simulation, base Stata. Model y = 1 + 0.5x + 0.3D + 0.4(x·D) + e with x ~ N(3, 1), D ~ Bernoulli(0.5), N = 500, 500 replications. Records the interaction t-statistic and VIF with raw and demeaned x. Runtime under one minute. Run from this folder: `do callen_centering.do`.
- `callen_mediator.do` — pure simulation, base Stata. Structural model rating = −0.6 earn + 0.5u + v, premium = −0.4 earn + rating + k·u + w, with k = 0 (design A) and k = 0.5 (design B); N = 2,000, 500 replications. Total effect −1.0, direct effect −0.4. Runtime under one minute. Run from this folder: `do callen_mediator.do`.
- `figures/fig1_centering.png` — interaction t-statistics with raw versus demeaned x (identical to 7 decimals; every point on the 45-degree line) and the interaction VIF (mean 11.1 raw versus 2.0 demeaned).
- `figures/fig2_mediator.png` — mean estimated earnings coefficient: −1.00 without the rating control (both designs), −0.40 with the control in design A, −0.28 with the control in design B (analytical value −0.280).

## What to look for in the output

- `callen_centering.do` prints the maximum absolute difference in the interaction t-statistic and standard error across replications (about 5e-7 and 8e-9, floating-point noise) alongside the mean VIFs. The VIF is a property of the columns of the design matrix; the standard error is a property of the parameter.
- `callen_mediator.do` prints the four means and the analytical population coefficient for design B: bias = 0.5 × cov(earn, u | rating) / var(earn | rating) = 0.5 × (0.3/1.61) / (1 − 0.36/1.61) = 0.120, so the coefficient is −0.280.

## Key readings

| Role | Paper |
|---|---|
| The focal paper | Callen, Livnat, and Segal (2009, *TAR*) |
| Centering and interactions | Aiken and West (1991); Echambadi and Hess (2007, *Marketing Science*); Kromrey and Foster-Johnson (1998, *Educational and Psychological Measurement*) |
| Bad controls and mediators | Rosenbaum (1984, *JRSS A*); Angrist and Pischke (2009, ch. 3); Cinelli, Forney, and Pearl (2022, *Sociological Methods & Research*) |
| The two-stage residual design | Beaver, Griffin, and Landsman (1982, *JAE*); Christie, Kennelley, King, and Schaefer (1984, *JAE*); Jennings (1990, *TAR*) |
