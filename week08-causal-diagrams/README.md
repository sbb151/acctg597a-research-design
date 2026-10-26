# Week 8: Controls and Designs as Causal Claims

**Paper:** Gow, I. D., D. F. Larcker, and P. C. Reiss (2016), "Causal Inference in Accounting Research," *Journal of Accounting Research* 54(2): 477-523.

## The design problem

Gow, Larcker, and Reiss coded every 2014 paper in the *Journal of Accounting Research*, *The Accounting Review*, and the *Journal of Accounting and Economics*: of 106 papers using observational data, 91 drew a causal inference, and 14 of those 91 relied on quasi-experimental methods (pp. 480-481, 485). The rest rely on control variables, and "the credibility of these assumptions is rarely explicitly addressed" (p. 481). Their Figure 1 (p. 483) makes the assumption visible: the same control Z is a confounder (Z causes X and Y), a mediator (X causes Z, which causes Y), or a collider (X and Y both cause Z) depending only on the direction of two arrows. Conditioning on Z is right in the first case and wrong in the other two, and "the inclusion of 'controls' that are mediators or colliders will generally lead to bias" (p. 484). Their running example is Larcker, Richardson, and Tuna (2007): governance factors, firm outcomes, and a list of controls whose roles are not stated; leverage, which governance and performance both plausibly affect, is their example of a collider (p. 485).

Quasi-experimental designs restate the same kind of claim as an absence. An instrument is valid only if no path runs from it to the outcome except through the treatment (GLR Figure 2, on Armstrong, Gow, and Larcker 2013), and "no statistical test allows the researcher to verify that their instruments satisfy the exclusion restriction" (p. 489). Difference-in-differences recovers a causal effect only if nothing that selected the treated group also moves their outcome's change; "if treatment assignment is driven by unobserved confounding variables, then DD and fixed-effect estimates of causal effects will be biased and inconsistent" (p. 493). The modules below make each claim concrete with a known truth.

## The simulation modules

`causal_diagrams.do` (base Stata, about a minute; the Monte Carlo in Module 3 dominates) runs three self-contained modules. In every design the causal effect of X on Y is 0.30 by construction. Identity checks are printed in the log and should be zero to eight decimals.

1. **One third variable, three roles.** Confounder: Z ~ N(0,1), X = 0.6 Z + u, Y = 0.30 X + 0.50 Z + e; the slope of Y on X is about 0.60 and falls to about 0.31 once Z is held fixed (omitted-variable identity exact). Mediator: M = 0.5 X + u, Y = 0.10 X + 0.40 M + e; the slope of Y on X is the total effect, about 0.31, and holding M fixed reports the direct effect, about 0.12 (c = c' + ab exact). Collider: Y = 0.30 X + e, C = 0.7 X + 0.7 Y + v; the slope of Y on X is about 0.31 without C and about -0.18 with C, against a closed form of -0.201. Output: `tables/tab_third_variable.tex`.
2. **Selection on a collider.** Keeping the observations with the highest values of C, from 100 to 10 percent of the sample, bends the retained-sample slope from about 0.31 through zero to about -0.06 in the top decile. Output: `figures/fig1_collider_selection.png`.
3. **A small exclusion violation and a weak instrument.** X = pi Z + 0.5 W + v and Y = 0.30 X + 0.05 Z + 0.5 W + e with W unobserved, so OLS is confounded (about 0.43 to 0.50) and the instrument has a direct effect on Y of one sixth of beta. Over 1,000 draws of N = 1,000 per point, the median IV estimate is about 0.35 at pi = 1 (first-stage F about 800) and about 0.78 at pi = 0.1 (F about 9), further from the truth than OLS, because the population IV slope is beta + gamma/pi. The 2SLS slope equals the reduced-form slope over the first-stage slope (identity exact). Output: `figures/fig2_iv_weak_instrument.png`.

The numbers above are from the Stata run (seed 2016) whose log, `causal_diagrams.log`, is included. The Week 8 deck was built while the Stata license on the build machine had lapsed, so its captions quote a NumPy replication of the same data-generating processes; the two differ in the second decimal for the collider slope and the weak-instrument points.

## Key readings

- Gow, Larcker, and Reiss (2016, *JAR*): the focal paper; the 2014 census, Pearl's causal diagrams as the organizing device, the audit of quasi-experimental methods, and the case for mechanisms, descriptive research, and structural models.
- Pearl (2009, *Causality*, 2nd ed.), as GLR's source: d-separation and the back-door criterion (GLR's appendix, pp. 517-519, states both).
- Larcker, Richardson, and Tuna (2007, *TAR*): GLR's running example of controls whose roles are unstated; their with-and-without-controls comparison and the "conservative" claim (pp. 484-485).
- Larcker and Rusticus (2010, *JAE*), as cited by GLR: instrument choice is a theoretical, not a statistical, exercise; GLR's 2014 IV papers "have paid little heed" (p. 491).
- Armstrong, Gow, and Larcker (2013, *JAR*) and Gow et al. (2013): the instrument GLR draw in Figure 2 and the finding that makes its exclusion restriction doubtful.
