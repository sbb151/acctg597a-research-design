********************************************************************************
* ACCTG 597A -- Research Design Module, Week 1
* The Reverse Regression Fallacy: Beaver, Lambert, and Ryan (1987, JAE)
*
* Purpose: Demonstrate by simulation that
*   (1) the reverse regression coefficient recovers the forward (structural)
*       coefficient only when R-squared = 1 (BLR's own Property 3: gamma*delta
*       = R^2), so with the R^2 ~ 0.13 typical of return-earnings regressions
*       the implied coefficient is inflated by a factor of 1/R^2 ~ 7.7;
*   (2) the "measurement error" BLR invoke lives in a latent construct
*       (permanent/economic income), not in observed earnings, so two
*       observationally equivalent worlds -- one with transitory earnings,
*       one with none -- produce identical regression output: delta is not
*       identified as a structural parameter;
*   (3) grouping data by the DEPENDENT variable mechanically inflates the
*       slope toward gamma/R^2 even with NO measurement error anywhere
*       (replicating BLR's Table 1 escalation 0.31 -> 2.30 in a world where
*       the true coefficient is 0.31);
*   (4) the same reverse-regression logic, inherited by Basu (1997), produces
*       spurious "asymmetric timeliness" in a world with zero conservatism.
*
* Requirements: base Stata (no user-written packages). Runtime ~1-2 minutes.
* All data are simulated; every true parameter is known by construction.
*
* Calibration: BLR (1987) Tables 1-2 -- forward gamma-hat = 0.31, mean
* R^2 = 0.13, reverse delta-hat = 0.59, grouped-by-G coefficient at the
* 5-portfolio level = 2.30, predicted asymptote gamma/R^2 = 2.38.
********************************************************************************

version 17
clear all
set more off
set seed 597

* All output lands next to this file. cd to the folder containing this .do
* file before running, or run via:  do reverse_regression.do
capture mkdir figures
capture log close
log using reverse_regression.log, replace text

********************************************************************************
* MODULE 1: The reverse regression recovers gamma only when R^2 = 1
*
* True model:  G = gamma*g + u,  gamma = 0.31 (BLR's own forward estimate).
* We vary var(u) to trace out R^2 from 0.05 to 0.99 and record:
*   - the forward estimate  gamma-hat        (unbiased at every R^2)
*   - the reverse estimate  delta-hat
*   - the implied forward coefficient 1/delta-hat
* BLR's Property (3), gamma*delta = R^2, implies 1/delta = gamma/R^2:
* the implied coefficient is inflated by exactly 1/R^2.
********************************************************************************

local gamma  = 0.31
local nobs   = 200000

tempname mem1
postfile `mem1' r2_target r2_hat gamma_fwd delta_rev implied theory ///
    using module1_results, replace

foreach r2 in 0.05 0.08 0.10 0.13 0.15 0.20 0.25 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 0.99 {
    quietly {
        clear
        set obs `nobs'
        * var(g) = 1.  R^2 = gamma^2 / (gamma^2 + var(u))  =>  solve for var(u)
        local varu = (`gamma'^2) * (1 - `r2') / `r2'
        gen g = rnormal(0, 1)
        gen G = `gamma'*g + rnormal(0, sqrt(`varu'))

        reg G g
        local gf  = _b[g]
        local r2h = e(r2)
        reg g G
        local d   = _b[G]

        post `mem1' (`r2') (`r2h') (`gf') (`d') (1/`d') (`gamma'/`r2')
    }
}
postclose `mem1'

use module1_results, clear
format gamma_fwd delta_rev implied theory %9.3f

di _n "MODULE 1: forward vs. implied-from-reverse estimates of gamma = `gamma'"
list r2_target gamma_fwd delta_rev implied theory, noobs sep(0)

* Figure 1: the implied coefficient explodes as R^2 falls
twoway (line theory r2_target, lcolor(gs10) lwidth(thick))                   ///
       (scatter implied r2_target, mcolor(maroon) msymbol(O))                ///
       (scatter gamma_fwd r2_target, mcolor(navy) msymbol(D)),               ///
    yline(0.31, lpattern(dash) lcolor(navy))                                 ///
    xline(0.13, lpattern(dot) lcolor(maroon))                                ///
    text(2.7 0.16 "BLR's world: R{sup:2} = 0.13,", place(e) size(small))      ///
    text(2.45 0.16 "implied coefficient = 0.31/0.13 = 2.38", place(e) size(small)) ///
    ytitle("Estimate of {&gamma}")                                           ///
    ylabel(0(1)6, angle(horizontal) format(%2.1f))                           ///
    xtitle("R{sup:2} of the return-earnings relation")                       ///
    title("The reverse regression is unbiased only when R{sup:2} = 1")       ///
    subtitle("True {&gamma} = 0.31; implied coefficient 1/{&delta} = {&gamma}/R{sup:2}") ///
    legend(order(2 "Implied from reverse regression (1/{&delta})"            ///
                 3 "Forward OLS ({&gamma}-hat)"                              ///
                 1 "Theory: {&gamma}/R{sup:2}") ring(0) pos(1) cols(1))      ///
    scheme(s1color)
graph export figures/fig1_implied_vs_r2.png, replace width(2000)

********************************************************************************
* MODULE 2: Two worlds, one regression -- delta is not identified
*
* World A (BLR taken literally): observed earnings changes contain a large
*   transitory component; prices respond only to the permanent part.
*       g = perm + trans,   G = perm + noise(small)
* World B (no transitory earnings at all): earnings are what they are --
*   the market reacts to all of g with ERC gamma0 -- but returns also carry
*   abundant NON-earnings information (discount rates, growth options, other
*   news). "Measurement error" exists only in the researcher's latent
*   construct (economic income), not in observed earnings.
*       g = perm,           G = gamma0*perm + noise(large)
*
* Both worlds are calibrated to BLR's moments: gamma-hat = 0.31, R^2 = 0.13.
* The forward, reverse, and grouped regressions are IDENTICAL across worlds.
* BLR interpret delta-hat as 1/(1-theta), the transitory-earnings parameter.
* World B has theta = 0 by construction and produces the same delta-hat.
********************************************************************************

* ---- World A: transitory earnings, clean returns --------------------------
* g = perm + trans with var(perm) = 0.31, var(trans) = 0.69 (so var(g) = 1),
* G = perm + e_A. Then cov(G,g) = var(perm) = 0.31 = gamma-hat. Choose
* var(e_A) so that R^2 = cov^2/(var(g)*var(G)) = 0.13 => var(G) = 0.739.
quietly {
    clear
    set obs `nobs'
    gen perm  = rnormal(0, sqrt(0.31))
    gen trans = rnormal(0, sqrt(0.69))
    gen g     = perm + trans
    gen G     = perm + rnormal(0, sqrt(0.739 - 0.31))
}
di _n "MODULE 2, World A: transitory earnings (theta > 0), little other news"
reg G g
local gfA = _b[g]
local r2A = e(r2)
reg g G
local dA  = _b[G]

* ---- World B: NO transitory earnings, noisy returns -----------------------
* g = perm (all of it priced, theta = 0), G = 0.31*g + e_B with var(e_B)
* chosen so R^2 = 0.13: var(e_B) = gamma^2*(1-R^2)/R^2 = 0.643.
quietly {
    clear
    set obs `nobs'
    gen g = rnormal(0, 1)
    gen G = 0.31*g + rnormal(0, sqrt(0.643))
}
di _n "MODULE 2, World B: zero transitory earnings (theta = 0), abundant non-earnings news"
reg G g
local gfB = _b[g]
local r2B = e(r2)
reg g G
local dB  = _b[G]

di _n "MODULE 2 summary: the same regressions in two different economies"
di    "                       World A (theta>0)   World B (theta=0)"
di    "forward gamma-hat:        " %6.3f `gfA' "             " %6.3f `gfB'
di    "R-squared:                " %6.3f `r2A' "             " %6.3f `r2B'
di    "reverse delta-hat:        " %6.3f `dA'  "             " %6.3f `dB'
di    "BLR-implied theta:        " %6.3f 1-1/`dA' "             " %6.3f 1-1/`dB'
di _n "delta-hat cannot distinguish transitory earnings from non-earnings"
di    "return variation: it is not identified as a structural parameter."

********************************************************************************
* MODULE 3: Grouping by the dependent variable is mechanical inflation
*
* Use World B: G = 0.31*g + u with NO measurement error in g anywhere.
* Group observations into portfolios by G (the dependent variable), take
* within-portfolio MEDIANS (the BLM/BLR convention), and rerun the forward
* regression at each grouping level: 100, 50, 25, 10, 5 portfolios.
* BLR's own appendix result: gamma-hat_p -> gamma/R^2 = 0.31/0.13 = 2.38.
* Their Table 1 reports exactly this escalation (0.31 -> 2.30) and BLM read
* it as evidence about measurement error. Here there is NO measurement
* error -- the escalation is pure statistical mechanics: portfolios with
* high G are portfolios whose disturbances u happened to be high, so
* grouping on G sorts on the error term itself.
* Contrast: grouping by g (the regressor) leaves the slope untouched.
********************************************************************************

* Mimic BLR's actual design: 19 annual cross-sections of ~700 firms each
* (their annual N ranges 439-771); run the grouped regression year by year
* and average the slopes across years, exactly as BLM/BLR do. With ~700
* firms, coarse portfolios (140 obs each) diversify the disturbance almost
* fully -- the slope approaches gamma/R^2 -- while fine portfolios (7 obs
* each) diversify it only partially, which is what generates BLR's
* signature escalation pattern across grouping levels.
local nyears = 19
local ncs3   = 700

tempname mem3
postfile `mem3' year nports str12 groupvar slope r2 using module3_results, replace

forvalues y = 1/`nyears' {
    quietly {
        clear
        set obs `ncs3'
        gen g = rnormal(0, 1)
        gen G = 0.31*g + rnormal(0, sqrt(0.643))
        save module3_year, replace

        reg G g
    }
    post `mem3' (`y') (`ncs3') ("none") (_b[g]) (e(r2))

    foreach gvar in G g {
        foreach p in 100 50 25 10 5 {
            quietly {
                use module3_year, clear
                xtile port = `gvar', nq(`p')
                collapse (median) G g, by(port)
                reg G g
            }
            post `mem3' (`y') (`p') ("`gvar'") (_b[g]) (e(r2))
        }
    }
}
postclose `mem3'

use module3_results, clear
collapse (mean) slope r2, by(groupvar nports)
gsort groupvar -nports
format slope r2 %9.3f
di _n "MODULE 3: mean forward slope across `nyears' simulated years of N = `ncs3'"
di    "          (true gamma = 0.31 in every year; NO measurement error anywhere)"
di    "          BLR Table 1, grouped by G: 0.31 -> 1.13 -> 1.58 -> 1.83 -> 2.07 -> 2.30"
list groupvar nports slope r2, noobs sepby(groupvar)
save module3_results, replace

* Figure 2: the escalation replicated without any measurement error,
* with BLR's actual Table 1 estimates overlaid for comparison
preserve
keep if groupvar != "none"
gen x = .
replace x = 1 if nports == 100
replace x = 2 if nports == 50
replace x = 3 if nports == 25
replace x = 4 if nports == 10
replace x = 5 if nports == 5
* BLR (1987), Table 1, grouped by G: reported mean coefficients
gen blr = .
replace blr = 1.13 if x == 1 & groupvar == "G"
replace blr = 1.58 if x == 2 & groupvar == "G"
replace blr = 1.83 if x == 3 & groupvar == "G"
replace blr = 2.07 if x == 4 & groupvar == "G"
replace blr = 2.30 if x == 5 & groupvar == "G"
twoway (connected slope x if groupvar == "G", mcolor(maroon) lcolor(maroon))  ///
       (connected blr x if groupvar == "G",                                   ///
            mcolor(maroon) msymbol(Th) lcolor(maroon) lpattern(shortdash))    ///
       (connected slope x if groupvar == "g", mcolor(navy) lcolor(navy)),     ///
    yline(0.31, lpattern(dash) lcolor(navy))                                  ///
    yline(2.38, lpattern(dash) lcolor(maroon))                                ///
    text(0.45 4.4 "true {&gamma} = 0.31", color(navy) size(small))            ///
    text(2.47 4.2 "asymptote {&gamma}/R{sup:2} = 2.38", color(maroon) size(small)) ///
    xlabel(1 "100" 2 "50" 3 "25" 4 "10" 5 "5")                                ///
    xtitle("Number of portfolios (coarser grouping {&rarr})")                 ///
    ytitle("Mean grouped-regression slope")                                   ///
    ylabel(0.5(0.5)2.5, angle(horizontal) format(%2.1f))                      ///
    title("Grouping by the dependent variable inflates the slope")            ///
    subtitle("Simulation has zero measurement error, yet reproduces BLR's escalation") ///
    legend(order(1 "Simulated: grouped by G (dependent variable)"             ///
                 2 "BLR Table 1, actual: grouped by G"                        ///
                 3 "Simulated: grouped by g (regressor)")                     ///
           ring(0) pos(11) cols(1) size(small))                               ///
    scheme(s1color)
graph export figures/fig2_grouping_escalation.png, replace width(2000)
restore

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
* Clean up intermediate datasets, keep results + figures + log
********************************************************************************
capture erase module3_year.dta

di _n "Done. Figures are in figures/; results datasets: module1_results.dta,"
di    "module3_results.dta, module4_results.dta; log: reverse_regression.log"

log close
