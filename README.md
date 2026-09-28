# ACCTG 597A: Research Design Modules

Simulation code accompanying the research-design modules in ACCTG 597A (PhD seminar, Penn State Smeal College of Business). Each week's folder demonstrates a research design problem that arose in the accounting literature and was institutionalized over time by citation to prior work that incorporated the flawed design.

Every `.do` file is **self-contained**: it generates its own simulated data, so the true parameter is known by construction and any bias is visible directly. No datasets, no WRDS access, and no user-written packages are required — base Stata (version 17 or later) is sufficient. Each folder's README explains the design flaw, the simulation modules, and the trail of econometric and accounting papers around it.

## Modules

| Week | Folder | Design flaw | Focal paper |
|---|---|---|---|
| 1 | `week01-reverse-regression/` | Reverse regression as a "fix" for errors-in-variables; grouping by the dependent variable | Beaver, Lambert, and Ryan (1987, *JAE*) |
| 1 | `week01-valuation-models/` | Certainty-world valuation models (MV = E/r, MV = BV) implicitly assumed by ERC and value-relevance regressions | Easton and Zmijewski (1989, *JAE*) |
| 3 | `week03-bgl-two-stage/` | The two-stage residual regression is one multiple regression read twice; orthogonalizing does not reduce collinearity; Simpson's paradox in pooled slopes | Beaver, Griffin, and Landsman (1982, *JAE*) |
| 3 | `week03-callen-controls/` | Demeaning does not reduce interaction-term multicollinearity; controlling for a mediator (the credit rating) changes the estimand | Callen, Livnat, and Segal (2009, *TAR*) |
| 3 | `week03-jennings-partial-derivatives/` | A coefficient is a partial derivative in a direction fixed by the other regressors; with linked regressors (an identity, clean surplus, a first stage, an interaction) the marginal effect is a linear combination of coefficients | Jennings (1990, *TAR*) |
| 4 | `week04-basu-asymmetric-timeliness/` | The asymmetric timeliness coefficient is nonzero under no conservatism: reversing the returns-earnings regression adds a variance-ratio bias and partitioning on the sign of returns adds a truncation bias, of unknown sign and size (Dietrich, Muller, and Riedl 2007) | Basu (1997, *JAE*) |
| 4 | `week04-bll-r2-incomparability/` | R-squared is not comparable across samples: it mixes the slope with the dispersion of the regressor, so scale controls do not restore comparability (Gu 2007) | Brown, Lo, and Lys (1999, *JAE*) |
| 4 | `week04-abk-heteroskedasticity/` | Deflation by a scale variable is weighted least squares under one assumed form of heteroskedasticity, correctly sized and efficient only when that form holds; White (1980) standard errors are consistent under any form, with a finite-sample leverage caveat that HC3 addresses | Aboody, Barth, and Kasznik (1999, *JAE*) |
| 5 | `week05-sloan-mishkin-test/` | The Mishkin test is, in large samples, a regression of abnormal returns on lagged information: its unconstrained estimates reproduce OLS, an omitted mispriced correlate shifts the starred coefficient, the implied misperception varies inversely with the earnings response coefficient, and pooled likelihood ratio statistics can overstate significance (Kraft, Leone, and Wasley 2007; Lewellen 2010) | Sloan (1996, *TAR*) |
| 5 | `week05-parallel-trends/` | Parallel trends restricts the untreated outcome of treated firms after treatment, which is never observed, so the assumption is not testable: a pre-trend test passes when the assumption fails (an observationally equivalent world), rejects when it holds (a transitory pre-period shock), and has low power against a drift (Roth 2022; Kahn-Lang and Lang 2020) | Schafhäutle and Veenman (2024, *TAR*) |
| 6 | `week06-mediation-moderation/` | The mediation decomposition c = c' + ab is an OLS identity, not a test: the indirect effect ab has a skewed sampling distribution (a percentile bootstrap interval over the Sobel test), a null total effect does not preclude it, the moderation coefficient b1 is a conditional effect at W = 0 and centering only relabels it, moderated mediation is tested through the index a3 b, and an unobserved variable in both M and Y or measurement error in M moves the split between direct and indirect paths while the total effect of an exogenous X stays unbiased (Hayes and Rockwood 2017; Lennox and Payne-Mann 2026) | Bhattacharya, Ecker, Olsson, and Schipper (2012, *TAR*) |

## Lecture decks

Slides for the in-class research-design modules are posted under `decks/` as PDFs.

| Week | Deck | Module |
|---|---|---|
| 1 | [`decks/week01/week1_reverse_regression.pdf`](decks/week01/week1_reverse_regression.pdf) | The Reverse Regression Fallacy (Beaver, Lambert, and Ryan 1987) |
| 1 | [`decks/week01/week1_valuation_models.pdf`](decks/week01/week1_valuation_models.pdf) | Valuation Models in Disguise (Easton and Zmijewski 1989) |
| 1 | [`decks/week01/week1_kl_borrowing.pdf`](decks/week01/week1_kl_borrowing.pdf) | Borrowed Technology (Kormendi and Lipe 1987) |
| 1 | [`decks/week01/week1_bmw_fact.pdf`](decks/week01/week1_bmw_fact.pdf) | A Fact Worth Explaining (Beaver, McNichols, and Wang 2020) |
| 2 | [`decks/week02/week2_ohlson_model.pdf`](decks/week02/week2_ohlson_model.pdf) | The Theory the Regressions Were Missing (Ohlson 1995) |
| 2 | [`decks/week02/week2_dechow_accruals.pdf`](decks/week02/week2_dechow_accruals.pdf) | The Value Added by Accountants (Dechow 1994) |
| 2 | [`decks/week02/week2_vuolteenaho_returns.pdf`](decks/week02/week2_vuolteenaho_returns.pdf) | What Drives Firm-Level Stock Returns? (Vuolteenaho 2002) |
| 2 | [`decks/week02/week2_michaely_signaling.pdf`](decks/week02/week2_michaely_signaling.pdf) | Signaling Safety (Michaely, Rossi, and Weber 2021) |
| 3 | [`decks/week03/week3_bgl_two_stage.pdf`](decks/week03/week3_bgl_two_stage.pdf) | Two Stages, One Regression (Beaver, Griffin, and Landsman 1982) |
| 3 | [`decks/week03/week3_callen_cds.pdf`](decks/week03/week3_callen_cds.pdf) | Centering, Mediators, and Residuals (Callen, Livnat, and Segal 2009) |
| 3 | [`decks/week03/week3_jennings_partials.pdf`](decks/week03/week3_jennings_partials.pdf) | Coefficients as Partial Derivatives (Jennings 1990) |
| 4 | [`decks/week04/week4_basu_asymmetric_timeliness.pdf`](decks/week04/week4_basu_asymmetric_timeliness.pdf) | Reversal and Truncation (Basu 1997) |
| 4 | [`decks/week04/week4_bll_r2_incomparability.pdf`](decks/week04/week4_bll_r2_incomparability.pdf) | R-Squared as a Sample Statistic (Brown, Lo, and Lys 1999) |
| 4 | [`decks/week04/week4_abk_heteroskedasticity.pdf`](decks/week04/week4_abk_heteroskedasticity.pdf) | Heteroskedasticity of Known and Unknown Form (Aboody, Barth, and Kasznik 1999) |
| 5 | [`decks/week05/week5_sloan_mishkin_test.pdf`](decks/week05/week5_sloan_mishkin_test.pdf) | The Mishkin Test as a Return Regression (Sloan 1996) |
| 5 | [`decks/week05/week5_sv_parallel_trends.pdf`](decks/week05/week5_sv_parallel_trends.pdf) | Parallel Trends as a Counterfactual Assumption (Schafhäutle and Veenman 2024) |
| 6 | [`decks/week06/week6_mediation_moderation.pdf`](decks/week06/week6_mediation_moderation.pdf) | Direct, Indirect, and Conditional Effects (Bhattacharya, Ecker, Olsson, and Schipper 2012) |

## Running the code

```bash
cd week01-reverse-regression
stata-mp -e do reverse_regression.do
```

or open the `.do` file in the Do-file Editor and run it. Figures are written to each folder's `figures/` subdirectory; a full log is saved alongside the code.

## Instructor

Sam Bonsall, Penn State Smeal College of Business
