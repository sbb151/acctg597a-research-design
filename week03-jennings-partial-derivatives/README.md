# Week 3: Coefficients as Partial Derivatives

**Paper:** Jennings (1990), "A Note on Interpreting 'Incremental Information Content'," *The Accounting Review* 65(4): 925-932.

## The design problem

Mid-1980s association studies attached the label "incremental information content" to any nonzero coefficient in a regression of long-window returns on the unexpected portions of two or more income variables. Jennings shows that the label conceals two hypotheses. A zero restriction on a component coefficient (CFO given TA) tests the *composition* of income: does the component add to the informativeness of earnings? An equality restriction between component coefficients, which is the same number as the coefficient on CFO given INC, tests the *disclosure* of income: are the components valued differently, so that separate reporting is useful? Because INC = CFO + TA, the two specifications are linear reparameterizations of each other, and "of the six coefficients that result from the three alternative regression specifications, only three are unique" (p. 927).

The general principle is that a regression coefficient is a partial derivative taken while the other regressors are held fixed. When the regressors are linked by an identity, a clean-surplus relation, a first-stage regression, or an interaction term, the marginal effect that answers the accounting question is a linear combination of coefficients, and the hypothesis test belongs on that combination (`lincom`, `test`), not on a single coefficient read off the table.

## The simulation modules

`jennings_partial_derivatives.do` (base Stata, a few seconds) runs four self-contained modules. Identity checks are printed in the log and should be zero to six decimals.

1. **Jennings' identity.** Two economies with 1,000 firms: equal valuation (RET = 0.30 CFO + 0.30 TA + e) and differential valuation (0.45 CFO + 0.15 TA). Each is estimated as eq. (2), RET on CFO and TA, and eq. (3), RET on CFO and INC. In every case gamma1 = theta1 - theta2 and gamma2 = theta2 exactly, with identical residuals and R-squared, and the F test of theta1 = theta2 equals the F test of gamma1 = 0. In the equal-valuation economy CFO is strongly informative (t about 6.7 in eq. 2) while its eq. (3) coefficient is about zero. Output: `tables/tab_reparam.tex`.
2. **Interaction terms.** RET = b0 + b1 E + b2 Size + b3 E x Size + e with Size roughly 2 to 10. The uncentered b1 is the earnings response at Size = 0, an extrapolation outside the data (about -0.36, t about -1.1); after mean centering, b1 is the response at the sample mean (about 1.58, t about 20.7). The figure plots b1 + b3 Size with a 95 percent band built from e(V). Output: `figures/fig1_interaction.png`.
3. **Ohlson levels regression under clean surplus.** MV = 2 + 5.0 E + 0.8 BV_t + e with BV_t = BV_t-1 + E - Div. The coefficient on E holding BV_t fixed is alpha1 (a dollar of earnings matched by a dollar of dividends); a retained dollar moves value by alpha1 + alpha2. Regressing MV on E and K = BV_t-1 - Div is an exact reparameterization and returns alpha1 + alpha2 on E with the `lincom` t-statistic. (A three-regressor version with BV_t-1 and Div entered separately is not exact, because it does not impose equal and opposite coefficients.) Output: `tables/tab_ohlson.tex`.
4. **BGL (1982) revisited.** On the seed-228 data of `week03-bgl-two-stage`, the derivative of returns with respect to HC is b1 with PRE held fixed (the Panel (B) residual coefficient), b1 + b2 with the depreciation adjustment D = PRE - HC held fixed, and b1 + c b2 with the first-stage residual Z held fixed (the Panel (A) coefficient on HC, which is also the bivariate slope). Three derivatives, one regression, each with a `lincom` t-statistic. Output: `tables/tab_bgl.tex`.

## Key readings

- Jennings (1990, *TAR*): the focal note; composition versus disclosure, the identity substitution, and the re-tests of Rayburn (1986) and Bowen, Burgstahler, and Daley (1987).
- Christie, Kennelley, King, and Schaefer (1984, *JAE*): the same reparameterization result, used to show that the BGL two-stage design is one multiple regression; cited in Jennings' footnote 4.
- Beaver, Griffin, and Landsman (1982, *JAE*): the two-stage residual regression revisited in Module 4 (see `week03-bgl-two-stage/`).
- Ohlson (1995, *CAR*): the clean-surplus relation (A2a) behind Module 3.
- Rayburn (1986, *JAR*) and Bowen, Burgstahler, and Daley (1987, *TAR*): the studies whose published tables Jennings re-reads.
- Biddle, Seow, and Siegel (1995, *CAR*): the later formalization of relative versus incremental information content.
