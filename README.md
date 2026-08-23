# ACCTG 597A: Research Design Modules

Simulation code accompanying the research-design modules in ACCTG 597A (PhD seminar, Penn State Smeal College of Business). Each week's folder demonstrates a research design problem that arose in the accounting literature and was institutionalized over time by citation to prior work that incorporated the flawed design.

Every `.do` file is **self-contained**: it generates its own simulated data, so the true parameter is known by construction and any bias is visible directly. No datasets, no WRDS access, and no user-written packages are required — base Stata (version 17 or later) is sufficient. Each folder's README explains the design flaw, the simulation modules, and the trail of econometric and accounting papers around it.

## Modules

| Week | Folder | Design flaw | Focal paper |
|---|---|---|---|
| 1 | `week01-reverse-regression/` | Reverse regression as a "fix" for errors-in-variables; grouping by the dependent variable | Beaver, Lambert, and Ryan (1987, *JAE*) |
| 4 | `week04-basu-asymmetric-timeliness/` | Spurious asymmetric timeliness from the piecewise reverse regression | Basu (1997, *JAE*) |

## Running the code

```bash
cd week01-reverse-regression
stata-mp -e do reverse_regression.do
```

or open the `.do` file in the Do-file Editor and run it. Figures are written to each folder's `figures/` subdirectory; a full log is saved alongside the code.

## Instructor

Sam Bonsall, Penn State Smeal College of Business
