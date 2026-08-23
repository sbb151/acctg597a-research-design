********************************************************************************
* ACCTG 597A -- Research Design Module, Week 4
* The Basu (1997) Legacy: Spurious Asymmetric Timeliness from the
* Piecewise Reverse Regression
*
* Companion to the Week 1 module (week01-reverse-regression/), which shows
* that the reverse regression of earnings on returns -- introduced by
* Beaver, Lambert, and Ryan (1987, JAE) -- recovers the forward coefficient
* only when R-squared = 1 and that grouping by the dependent variable is
* mechanically biased. Basu (1997, JAE) builds the conditional-conservatism
* literature on the same earnings-on-returns direction, adding a piecewise
* split on the sign of returns. This module shows that design produces
* spurious "asymmetric timeliness" in a world with zero conservatism.
*
* Requirements: base Stata (no user-written packages). Runtime ~1 minute.
* All data are simulated; every true parameter is known by construction.
********************************************************************************

version 17
clear all
set more off
set seed 597

* All output lands next to this file. cd to the folder containing this .do
* file before running, or run via:  do basu_spurious_timeliness.do
capture mkdir figures
capture log close
log using basu_spurious_timeliness.log, replace text

********************************************************************************
* MODULE 4: The Basu (1997) legacy -- spurious asymmetric timeliness
*
* Basu's conservatism regression is a piecewise REVERSE regression:
*     g = b0 + b1*D + b2*G + b3*(D*G) + e,   D = 1(G < 0),
* and b3 > 0 is read as conservatism (earnings timelier for bad news).
*
* Simulate a world with ZERO conservatism: earnings respond to news
* symmetrically. But let the non-earnings component of returns be
* right-skewed (limited liability; growth-option and discount-rate news
* produce occasional large positive returns unrelated to current earnings).
* Then negative returns are disproportionately NEWS-driven while large
* positive returns are disproportionately NOISE-driven, so the reverse
* slope cov(g,G)/var(G) is steeper on the negative side: b3 > 0 with no
* conservatism whatsoever. The Basu coefficient inherits every pathology
* of the reverse regression: it is a statement about the variance
* composition of returns, not about accounting.
*
* Placebo: with SYMMETRIC non-earnings noise, b3 = 0 -- confirming the
* skewness of return noise, not accounting behavior, drives the result.
********************************************************************************

* Monte Carlo at a realistic cross-section size (BLR's annual N ~ 700)
local reps  = 1000
local ncs   = 700

tempname mem4
postfile `mem4' rep b3_skew t3_skew b3_symm t3_symm using module4_results, replace

forvalues r = 1/`reps' {
    quietly {
        clear
        set obs `ncs'
        gen news = rnormal(0, 1)
        * earnings: symmetric response to news + idiosyncratic noise
        gen g = 0.31*news + rnormal(0, 0.5)

        * returns, version 1: right-skewed non-earnings noise
        * (standardized chi-square(2): mean 0, variance 1, skewness 2)
        gen zskew = (rchi2(2) - 2)/2
        gen G1 = news + 1.4*zskew

        * returns, version 2: symmetric non-earnings noise, same variance
        gen G2 = news + 1.4*rnormal(0, 1)

        * Basu regression on each
        gen D1  = (G1 < 0)
        gen DG1 = D1*G1
        reg g D1 G1 DG1
        local b3s = _b[DG1]
        local t3s = _b[DG1]/_se[DG1]

        gen D2  = (G2 < 0)
        gen DG2 = D2*G2
        reg g D2 G2 DG2
        local b3n = _b[DG2]
        local t3n = _b[DG2]/_se[DG2]

        post `mem4' (`r') (`b3s') (`t3s') (`b3n') (`t3n')
    }
}
postclose `mem4'

use module4_results, clear
gen rej_skew = abs(t3_skew) > 1.96
gen rej_symm = abs(t3_symm) > 1.96

di _n "MODULE 4: Basu asymmetric-timeliness coefficient in a world with NO conservatism"
di    "          (true b3 = 0 by construction; `reps' Monte Carlo samples of N = `ncs')"
su b3_skew t3_skew rej_skew b3_symm t3_symm rej_symm

quietly su b3_skew
local mb3 : di %5.3f r(mean)
quietly su rej_skew
local rr  : di %4.1f 100*r(mean)

* Figure 3: sampling distribution of the spurious Basu coefficient
twoway (histogram b3_symm, width(0.008) color(navy%40))                       ///
       (histogram b3_skew, width(0.008) color(maroon%40)),                    ///
    xline(0, lpattern(dash) lcolor(black))                                    ///
    xtitle("Estimated asymmetric-timeliness coefficient (b{sub:3})")          ///
    ytitle("Density")                                                         ///
    ylabel(0(2)10, angle(horizontal) format(%2.0f))                           ///
    title("Spurious conservatism from skewed return noise")                   ///
    subtitle("True b{sub:3} = 0 in both worlds; mean spurious b{sub:3} = `mb3', rejection rate `rr'%") ///
    legend(order(2 "Right-skewed non-earnings return noise"                   ///
                 1 "Symmetric non-earnings return noise (placebo)")           ///
           ring(1) pos(6) cols(1) size(small))                                ///
    scheme(s1color)
graph export figures/fig3_spurious_basu.png, replace width(2000)

* Figure 4: one large sample, the piecewise fit that "finds" conservatism
quietly {
    clear
    set obs 50000
    gen news  = rnormal(0, 1)
    gen g     = 0.31*news + rnormal(0, 0.5)
    gen zskew = (rchi2(2) - 2)/2
    gen G     = news + 1.4*zskew
    gen D     = (G < 0)
    gen DG    = D*G
    reg g D G DG
    predict ghat
    * portfolio means for a readable scatter
    xtile bin = G, nq(50)
    collapse (mean) g G ghat, by(bin)
}
twoway (scatter g G, mcolor(gs8) msymbol(oh))                                 ///
       (line ghat G if G < 0, lcolor(maroon) lwidth(thick))                   ///
       (line ghat G if G >= 0, lcolor(navy) lwidth(thick)),                   ///
    xline(0, lpattern(dot) lcolor(black))                                     ///
    xtitle("Return (G)") ytitle("Earnings change (g)")                        ///
    title("The piecewise fit 'finds' conservatism")                           ///
    subtitle("True earnings timeliness is perfectly symmetric")               ///
    note("Portfolio means (50 bins) of 50,000 simulated observations")        ///
    ylabel(, angle(horizontal) format(%2.1f))                                 ///
    legend(order(2 "Fitted slope, bad news (steeper)"                         ///
                 3 "Fitted slope, good news") ring(0) pos(5) cols(1))         ///
    scheme(s1color) scale(1.15)
graph export figures/fig4_basu_piecewise.png, replace width(2000)

********************************************************************************
* Keep results + figures + log
********************************************************************************

di _n "Done. Figures are in figures/; results dataset: module4_results.dta;"
di    "log: basu_spurious_timeliness.log"

log close
