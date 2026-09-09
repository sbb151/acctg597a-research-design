********************************************************************************
* ACCTG 597A -- Week 3: Two Stages, One Regression
* Beaver, Griffin, and Landsman (1982, JAE), "The Incremental Information
* Content of Replacement Cost Earnings," and the two-stage residual regression.
*
* Self-contained simulation, base Stata only (version 17 or later). Run from
* this folder:  do bgl_two_stage.do
* Figures are written to figures/. Runtime about two minutes (the Monte Carlo
* in Module 3 dominates).
*
* Four modules:
*   1. Correlation tables hold nothing fixed. From BGL Table 3 (1977) the
*      partial correlation of PRE with return given HC is -0.02, although the
*      pairwise correlation is 0.30.
*   2. Both panels of BGL Table 4 are one regression read twice. On simulated
*      data calibrated to BGL, the Panel (A) residual coefficient equals the
*      one-regression coefficient on PRE, the Panel (B) residual coefficient
*      equals the one-regression coefficient on HC (same t-statistics), and
*      the un-orthogonalized coefficients are bivariate slopes.
*   3. Orthogonalizing does not reduce the standard error of the incremental
*      term. Monte Carlo across correlations: the two-stage and one-regression
*      sampling standard deviations coincide and follow the textbook formula.
*   4. Simpson's paradox. Identical year-by-year slopes can pool to a slope of
*      the opposite sign, and zero year slopes can pool to a positive slope,
*      so agreement between year and pooled columns is not evidence.
*
* Calibration: BGL Table 2 (sd of HC 0.216, sd of PRE 0.377, sd of return
* about 0.21), Table 3 (corr(HC, PRE) = 0.84, corr(HC, R) = 0.37,
* corr(PRE, R) = 0.30), Table 4 (one-regression coefficients 0.39 on HC and
* -0.02 on PRE, R-squared 0.14, 303 firms).
********************************************************************************

version 17
clear all
set more off
capture mkdir figures

********************************************************************************
* Module 1: pairwise versus partial correlations (recreated from Table 3)
********************************************************************************
scalar r_hr  = 0.37
scalar r_pr  = 0.30
scalar r_hp  = 0.84
scalar p_pre = (r_pr - r_hp*r_hr) / sqrt((1 - r_hp^2)*(1 - r_hr^2))
scalar p_hc  = (r_hr - r_hp*r_pr) / sqrt((1 - r_hp^2)*(1 - r_pr^2))
di _n "Module 1: partial corr(PRE, R | HC) = " %6.3f p_pre ///
      ";  partial corr(HC, R | PRE) = " %6.3f p_hc

clear
set obs 4
gen x     = .
gen value = .
gen kind  = .
replace x = 1 in 1
replace x = 2 in 2
replace x = 4 in 3
replace x = 5 in 4
replace value = r_pr  in 1
replace value = p_pre in 2
replace value = r_hr  in 3
replace value = p_hc  in 4
replace kind = 1 in 1
replace kind = 2 in 2
replace kind = 1 in 3
replace kind = 2 in 4
local lab_pre_pair : di %4.2f r_pr
local lab_pre_part : di %5.2f p_pre
local lab_hc_pair  : di %4.2f r_hr
local lab_hc_part  : di %4.2f p_hc

twoway (bar value x if kind == 1, barwidth(0.8) color(navy))               ///
       (bar value x if kind == 2, barwidth(0.8) color(maroon)),            ///
    text(0.335 1 "`lab_pre_pair'", color(navy) size(medium))               ///
    text(0.03  2 "`lab_pre_part'", color(maroon) size(medium))             ///
    text(0.405 4 "`lab_hc_pair'", color(navy) size(medium))                ///
    text(0.26  5 "`lab_hc_part'", color(maroon) size(medium))              ///
    text(0.47 1.5 "Replacement cost (PRE)", color(black) size(medsmall))   ///
    text(0.47 4.5 "Historical cost (HC)", color(black) size(medsmall))     ///
    xlabel(1 "pairwise" 2 "partial" 4 "pairwise" 5 "partial", noticks)     ///
    xscale(range(0.3 5.7))                                                 ///
    ylabel(-0.1(0.1)0.5, angle(horizontal) format(%3.1f))                  ///
    yline(0, lcolor(gs8) lwidth(thin))                                     ///
    xtitle("") ytitle("Correlation with security return")                  ///
    title("Pairwise 0.30 becomes partial -0.02 given HC")                  ///
    subtitle("Recreated from BGL (1982) Table 3, 1977")                    ///
    legend(off) scheme(s1color)
graph export figures/fig1_pairwise_partial.png, replace width(2000)

********************************************************************************
* Module 2: both panels of Table 4 are one regression read twice
********************************************************************************
clear
set seed 228
set obs 303
gen hc  = 0.134 + 0.216*rnormal()
gen pre = 0.148 + 1.46*(hc - 0.134) + 0.206*rnormal()
gen r   = 0.018 + 0.39*(hc - 0.134) - 0.02*(pre - 0.148) + 0.19*rnormal()

di _n "Module 2: simulated correlations (compare BGL Table 3: 0.84, 0.37, 0.30)"
corr hc pre r

* Panel (A): orthogonalize PRE against HC
reg pre hc
scalar cA = _b[hc]
predict zA, resid
reg r hc zA
scalar A_hc = _b[hc]
scalar A_z  = _b[zA]
scalar A_zt = _b[zA]/_se[zA]

* Panel (B): orthogonalize HC against PRE
reg hc pre
scalar cB = _b[pre]
predict zB, resid
reg r pre zB
scalar B_pre = _b[pre]
scalar B_z   = _b[zB]
scalar B_zt  = _b[zB]/_se[zB]

* One regression
reg r hc pre
scalar O_hc   = _b[hc]
scalar O_hct  = _b[hc]/_se[hc]
scalar O_pre  = _b[pre]
scalar O_pret = _b[pre]/_se[pre]

* Bivariate slopes
reg r hc
scalar biv_hc = _b[hc]
reg r pre
scalar biv_pre = _b[pre]

di _n "Identity checks (zero to machine precision):"
di "  Panel (A) residual coef - one-regression PRE coef : " %10.6f (A_z - O_pre)
di "  Panel (A) residual t    - one-regression PRE t    : " %10.6f (A_zt - O_pret)
di "  Panel (B) residual coef - one-regression HC coef  : " %10.6f (B_z - O_hc)
di "  Panel (B) residual t    - one-regression HC t     : " %10.6f (B_zt - O_hct)
di "  Panel (A) HC coef - (b_HC + c_A b_PRE)             : " %10.6f (A_hc - (O_hc + cA*O_pre))
di "  Panel (A) HC coef - bivariate slope of R on HC     : " %10.6f (A_hc - biv_hc)
di "  Panel (B) PRE coef - (b_PRE + c_B b_HC)            : " %10.6f (B_pre - (O_pre + cB*O_hc))
di "  Panel (B) PRE coef - bivariate slope of R on PRE   : " %10.6f (B_pre - biv_pre)

* Added-variable plot (Frisch-Waugh-Lovell): the one-regression coefficient
* on PRE, drawn as a picture
reg r hc
predict r_res, resid
reg pre hc
predict pre_res, resid
reg r pre
local b_raw : di %4.2f _b[pre]
reg r_res pre_res
local b_part : di %4.2f _b[pre_res]
local t_part : di %3.1f (_b[pre_res]/_se[pre_res])

twoway (scatter r pre, mcolor(navy%45) msize(small))                       ///
       (lfit r pre, lcolor(navy) lwidth(thick)),                           ///
    text(0.66 -0.95 "slope `b_raw'", color(navy) size(medium) placement(e)) ///
    xtitle("PRE (percentage change)") ytitle("Return")                     ///
    xlabel(-1(0.5)1.5, format(%3.1f)) ylabel(-0.6(0.3)0.6, angle(horizontal) format(%3.1f)) ///
    yscale(range(-0.7 0.75))                                               ///
    title("Raw: R on PRE") legend(off) scheme(s1color)                     ///
    name(g_raw, replace) nodraw
twoway (scatter r_res pre_res, mcolor(maroon%45) msize(small))             ///
       (lfit r_res pre_res, lcolor(maroon) lwidth(thick)),                 ///
    text(0.66 -0.95 "slope `b_part' (t = `t_part')", color(maroon) size(medium) placement(e)) ///
    xtitle("PRE residual, given HC") ytitle("Return residual, given HC")   ///
    xlabel(-1(0.5)1, format(%3.1f)) ylabel(-0.6(0.3)0.6, angle(horizontal) format(%3.1f)) ///
    yscale(range(-0.7 0.75))                                               ///
    title("Holding HC fixed: added-variable plot") legend(off) scheme(s1color) ///
    name(g_avp, replace) nodraw
graph combine g_raw g_avp, rows(1) xsize(10) ysize(4.4)                     ///
    title("The PRE association vanishes once HC is held fixed")             ///
    subtitle("Simulated data calibrated to BGL (1982) Tables 2-4, N = 303") ///
    scheme(s1color)
graph export figures/fig2_added_variable.png, replace width(2000)

********************************************************************************
* Module 3: orthogonalizing does not shrink the standard error
********************************************************************************
set seed 19841201
local N      = 303
local reps   = 400
local sig_e  = 0.19
local sd_hc  = 0.216
local sd_pre = 0.377

tempname sims
tempfile simsfile repfile curvefile
postfile `sims' rho sd_one sd_two se_formula using `simsfile', replace
di _n "Module 3: Monte Carlo (400 replications per correlation)"
foreach rho in 0 0.3 0.5 0.7 0.84 0.9 0.95 0.98 {
    tempname acc
    postfile `acc' b_one b_two using `repfile', replace
    forvalues i = 1/`reps' {
        quietly {
            clear
            set obs `N'
            gen hc  = `sd_hc'*rnormal()
            gen pre = `sd_pre'*(`rho'*(hc/`sd_hc') + sqrt(1 - `rho'^2)*rnormal())
            gen y   = 0.39*hc + `sig_e'*rnormal()
            reg y hc pre
            scalar b1 = _b[pre]
            reg pre hc
            predict z, resid
            reg y hc z
            scalar b2 = _b[z]
            post `acc' (b1) (b2)
        }
    }
    postclose `acc'
    quietly use `repfile', clear
    quietly sum b_one
    scalar s1 = r(sd)
    quietly sum b_two
    scalar s2 = r(sd)
    scalar sf = `sig_e'/(sqrt(`N')*`sd_pre'*sqrt(1 - `rho'^2))
    di "  rho = " %4.2f `rho' "  sd(one regression) = " %6.4f s1 ///
       "  sd(two-stage) = " %6.4f s2 "  formula = " %6.4f sf
    post `sims' (`rho') (s1) (s2) (sf)
}
postclose `sims'

clear
set obs 200
gen rho_c = 0.985*(_n - 1)/199
gen se_c  = `sig_e'/(sqrt(`N')*`sd_pre'*sqrt(1 - rho_c^2))
save `curvefile', replace
use `simsfile', clear
append using `curvefile'

twoway (line se_c rho_c, lcolor(gs10) lwidth(medthick))                     ///
       (scatter sd_one rho, mcolor(navy) msymbol(circle) msize(large))       ///
       (scatter sd_two rho, mcolor(maroon) msymbol(square_hollow) msize(large) mlwidth(medthick)), ///
    xline(0.84, lcolor(gs8) lpattern(dash))                                  ///
    text(0.150 0.82 "BGL: corr = 0.84, VIF = 3.4", color(black) size(medium) placement(w)) ///
    text(0.095 0.02 "navy circles: one regression", color(navy) size(medium) placement(e)) ///
    text(0.080 0.02 "maroon squares: two-stage residual coefficient", color(maroon) size(medium) placement(e)) ///
    text(0.065 0.02 "gray curve: textbook formula", color(gs8) size(medium) placement(e)) ///
    xtitle("Correlation between HC and PRE")                                 ///
    ytitle("Sampling s.d. of the incremental PRE coefficient")               ///
    xlabel(0(0.2)1, format(%3.1f)) ylabel(0(0.03)0.18, angle(horizontal) format(%4.2f)) ///
    yscale(range(0 0.185))                                                   ///
    title("Orthogonalizing leaves the standard error unchanged")             ///
    subtitle("Monte Carlo, 400 cross-sections of 303 firms per point")       ///
    legend(off) scheme(s1color) xsize(10) ysize(5.2)
graph export figures/fig3_collinearity_se.png, replace width(2000)

********************************************************************************
* Module 4: Simpson's paradox
********************************************************************************
set seed 19771231
clear
set obs 303
gen firm = _n
gen u  = 0.216*rnormal()
gen e  = 0.19*rnormal()
tempfile draws
save `draws', replace

* Panel (a): identical within-year slopes (0.36), pooled slope negative
use `draws', clear
expand 2
bysort firm: gen year = cond(_n == 1, 1977, 1978)
gen hc = cond(year == 1977, 0.10, 0.55) + u
gen r  = cond(year == 1977, 0.16, -0.16) + 0.36*u + e
reg r hc if year == 1977
local a77 : di %4.2f _b[hc]
reg r hc if year == 1978
local a78 : di %4.2f _b[hc]
reg r hc
local apool : di %4.2f _b[hc]
reg r hc i.year
local afe : di %4.2f _b[hc]
di _n "Module 4, panel (a): 1977 = `a77', 1978 = `a78', pooled = `apool', pooled with year effects = `afe'"

twoway (scatter r hc if year == 1977, mcolor(navy%40) msize(small))        ///
       (scatter r hc if year == 1978, mcolor(maroon%40) msize(small))      ///
       (lfit r hc if year == 1977, lcolor(navy) lwidth(thick))             ///
       (lfit r hc if year == 1978, lcolor(maroon) lwidth(thick))           ///
       (lfit r hc, lcolor(black) lwidth(thick) lpattern(dash)),            ///
    text(0.74 -0.55 "1977: slope `a77'", color(navy) size(medsmall) placement(e))   ///
    text(0.74 0.66 "1978: slope `a78'", color(maroon) size(medsmall) placement(e))  ///
    text(-0.72 0.55 "pooled: slope `apool'", color(black) size(medsmall) placement(e)) ///
    xtitle("HC (percentage change)") ytitle("Return")                       ///
    xlabel(-0.5(0.5)1.5, format(%3.1f)) ylabel(-0.6(0.3)0.6, angle(horizontal) format(%3.1f)) ///
    yscale(range(-0.8 0.78))                                                ///
    title("(a) Identical year slopes, negative pooled slope")               ///
    legend(off) scheme(s1color) name(g_a, replace) nodraw

* Panel (b): zero within-year slopes, pooled slope positive
use `draws', clear
expand 2
bysort firm: gen year = cond(_n == 1, 1977, 1978)
gen hc = cond(year == 1977, 0.10, 0.55) + u
gen r  = cond(year == 1977, -0.15, 0.16) + 0.0*u + e
reg r hc if year == 1977
local b77 : di %4.2f _b[hc]
reg r hc if year == 1978
local b78 : di %4.2f _b[hc]
reg r hc
local bpool : di %4.2f _b[hc]
reg r hc i.year
local bfe : di %4.2f _b[hc]
di "Module 4, panel (b): 1977 = `b77', 1978 = `b78', pooled = `bpool', pooled with year effects = `bfe'"

twoway (scatter r hc if year == 1977, mcolor(navy%40) msize(small))        ///
       (scatter r hc if year == 1978, mcolor(maroon%40) msize(small))      ///
       (lfit r hc if year == 1977, lcolor(navy) lwidth(thick))             ///
       (lfit r hc if year == 1978, lcolor(maroon) lwidth(thick))           ///
       (lfit r hc, lcolor(black) lwidth(thick) lpattern(dash)),            ///
    text(-0.76 -0.55 "1977: slope `b77'", color(navy) size(medsmall) placement(e))  ///
    text(-0.76 0.55 "1978: slope `b78'", color(maroon) size(medsmall) placement(e)) ///
    text(0.68 -0.55 "pooled: slope `bpool'", color(black) size(medsmall) placement(e)) ///
    xtitle("HC (percentage change)") ytitle("Return")                       ///
    xlabel(-0.5(0.5)1.5, format(%3.1f)) ylabel(-0.6(0.3)0.6, angle(horizontal) format(%3.1f)) ///
    yscale(range(-0.8 0.78))                                                ///
    title("(b) Zero year slopes, positive pooled slope")                    ///
    legend(off) scheme(s1color) name(g_b, replace) nodraw

graph combine g_a g_b, rows(1) xsize(10) ysize(4.4)                          ///
    title("Agreement between year and pooled slopes is not evidence")        ///
    subtitle("Simulated cross-sections, 303 firms per year")                 ///
    scheme(s1color)
graph export figures/fig4_simpson.png, replace width(2000)

di _n "Done. Figures written to figures/."
