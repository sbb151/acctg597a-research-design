********************************************************************************
* ACCTG 597A -- Week 6: Direct, Indirect, and Conditional Effects
* Bhattacharya, Ecker, Olsson, and Schipper (2012), "Direct and Mediated
* Associations among Earnings Quality, Information Asymmetry, and the Cost of
* Equity," The Accounting Review 87(2): 449-482; Hayes and Rockwood (2017),
* "Regression-based statistical mediation and moderation analysis in clinical
* research," Behaviour Research and Therapy 98: 39-57; Lennox and Payne-Mann
* (2026), "An Explanation of Path Analysis and Recommendations for Best
* Practice," Contemporary Accounting Research 43(1): 290-313.
*
* Self-contained simulation, base Stata only (version 17 or later; Module 2
* uses Mata, which ships with Stata). Run from this folder:
*     do mediation_moderation.do
* Figures are written to figures/ and booktabs table fragments to tables/.
* Runtime about three minutes (Module 2 runs 5,000 samples with 1,000
* bootstrap resamples each).
*
* Notation follows Hayes and Rockwood (2017): M = i_M + a X, Y = i_Y + c' X +
* b M, Y = i_Y + c X; moderation Y = i + b1 X + b2 W + b3 XW; conditional
* process M = a1 X + a2 W + a3 XW, Y = c' X + b M.
*
* Modules. In each, the true parameter is known by construction, so the
* property under study is visible directly.
*   1. Identity. M = 0.5 X + e_M; Y = 0.2 X + 0.8 M + e_Y (Design A) or
*      Y = -0.4 X + 0.8 M + e_Y (Design B), N = 5,000, seed 2012. c = c' + ab
*      to machine precision in both designs; Design B has ab = 0.38 with c
*      indistinguishable from zero; sem reproduces the OLS paths.
*   2. Bootstrap. 5,000 samples of N = 100 with ab = 0.0625, seed 2017.
*      Skewness of ab-hat about 0.78; Sobel interval coverage about 0.92
*      against 0.94 for the percentile bootstrap; then the one-sample recipe
*      with Stata's bootstrap command.
*   3. Conditional effect. Y = 0.15 X + 0.30 W + 0.25 XW + e, N = 500, seed
*      2017. theta(W) = b1 + b3 W; centering relabels b1 and leaves b3
*      unchanged; Johnson-Neyman boundaries near W = -0.89 and -0.21; a
*      median split is not a test of b3.
*   4. Moderated mediation. M = 0.5 X + 0.2 W + 0.3 XW + e_M, Y = 0.2 X +
*      0.6 M + e_Y, N = 500, seed 2017. omega(W) = (a1 + a3 W) b; index
*      a3 b about 0.22 with a percentile bootstrap interval; reduced-form
*      identities exact.
*   5. Path bias. 500 samples of N = 2,000, seed 2012. A shared unobservable
*      in M and Y or measurement error in M moves the total between the
*      direct and indirect paths (0.07/0.53 and 0.40/0.20 against the true
*      0.20/0.40) while c' + ab = c holds in every sample.
********************************************************************************

version 17
clear all
set more off
capture log close
log using mediation_moderation.log, replace text
capture mkdir figures
capture mkdir tables

********************************************************************************
* MODULE 1: the decomposition c = c' + ab holds exactly; a zero total effect does not preclude an indirect effect
********************************************************************************
********************************************************************************
* ACCTG 597A -- Week 6 (mediation and moderation module)
* The total effect decomposes exactly into direct and indirect components:
* c = c' + ab when one OLS sample estimates M on X and Y on X and M.
*
* Pure simulation, N = 5,000. One draw of X and e_M is shared by two designs,
* so the a path is identical and only the Y equation differs.
*   M = 0.5 X + e_M                          (a = 0.5)
*   Design A: Y = 0.2 X + 0.8 M + e_Y        (c' = 0.2, b = 0.8, ab = 0.4, c = 0.6)
*   Design B: Y = -0.4 X + 0.8 M + e_Y       (c' = -0.4, b = 0.8, ab = 0.4, c = 0.0)
* All shocks standard normal. Design B has a zero total effect and a nonzero
* indirect effect, so a test of c cannot serve as a gate for a test of ab.
* Identity checked: c - (c' + ab) = 0 to machine precision in both designs;
* the path coefficients from Stata's sem command equal the OLS coefficients.
* Requirements: base Stata. Output: tables/tab_mediation_moderation_identity.tex
********************************************************************************

clear all
set seed 2012

set obs 5000
gen x  = rnormal()
gen em = rnormal()
gen ey = rnormal()
gen m  = 0.5*x + em
gen yA =  0.2*x + 0.8*m + ey
gen yB = -0.4*x + 0.8*m + ey

* ---- Eq. (2) of Hayes and Rockwood: M on X (the a path), common to both designs ----
reg m x
scalar a   = _b[x]
scalar at  = _b[x]/_se[x]
scalar ase = _se[x]

foreach d in A B {
    * Eq. (3): Y on X and M (c' and b)
    reg y`d' x m
    scalar cp`d'   = _b[x]
    scalar cp`d't  = _b[x]/_se[x]
    scalar b`d'    = _b[m]
    scalar b`d't   = _b[m]/_se[m]
    scalar bse`d'  = _se[m]
    * Eq. (1): Y on X (the total effect c)
    reg y`d' x
    scalar c`d'    = _b[x]
    scalar c`d't   = _b[x]/_se[x]
    scalar ab`d'   = a*b`d'
    scalar sum`d'  = cp`d' + ab`d'
    scalar sob`d'  = sqrt(b`d'^2*ase^2 + a^2*bse`d'^2)
    scalar abz`d'  = ab`d'/sob`d'
}

di _n "Design A: a = " %6.3f a " b = " %6.3f bA " c' = " %6.3f cpA ///
      " ab = " %6.3f abA " (Sobel z = " %4.1f abzA ") c = " %6.3f cA " (t = " %4.1f cAt ")"
di    "Design B: a = " %6.3f a " b = " %6.3f bB " c' = " %6.3f cpB ///
      " ab = " %6.3f abB " (Sobel z = " %4.1f abzB ") c = " %6.3f cB " (t = " %4.1f cBt ")"

di _n "Identity checks (should be zero to machine precision):"
foreach d in A B {
    di "  Design `d': c - (c' + ab) = " %12.8f (c`d' - (cp`d' + ab`d'))
    assert abs(c`d' - (cp`d' + ab`d')) < 1e-6
}

* ---- Path analysis by maximum likelihood (Stata's sem) gives the same paths ----
sem (m <- x) (yA <- x m)
scalar a_sem  = _b[m:x]
scalar b_sem  = _b[yA:m]
scalar cp_sem = _b[yA:x]
estat teffects
di _n "  sem a*b minus OLS ab (Design A)   = " %12.8f (a_sem*b_sem - abA)
di    "  sem c' minus OLS c' (Design A)    = " %12.8f (cp_sem - cpA)
assert abs(a_sem*b_sem - abA) < 1e-6
assert abs(cp_sem - cpA) < 1e-6

* ---- Numbers for the table ----
foreach s in a cpA bA cA abA sumA cpB bB cB abB sumB {
    local f`s' : di %4.2f `s'
    local f`s' = strtrim("`f`s''")
}
foreach s in at cpAt bAt cAt abzA cpBt bBt cBt abzB {
    local f`s' : di %3.1f `s'
    local f`s' = strtrim("`f`s''")
}

capture mkdir tables
file open tab using tables/tab_mediation_moderation_identity.tex, write replace
file write tab "\begin{tabular}{@{}lcc@{}}" _n
file write tab "\toprule" _n
file write tab " & Design A & Design B \\" _n
file write tab "\midrule" _n
file write tab "\$a\$: \$M\$ on \$X\$ & `fa' (`fat') & `fa' (`fat') \\" _n
file write tab "\$b\$: \$Y\$ on \$M\$, \$X\$ held fixed & `fbA' (`fbAt') & `fbB' (`fbBt') \\" _n
file write tab "\$c'\$: \$Y\$ on \$X\$, \$M\$ held fixed & `fcpA' (`fcpAt') & `fcpB' (`fcpBt') \\" _n
file write tab "\midrule" _n
file write tab "\$ab\$ (indirect) & `fabA' [`fabzA'] & `fabB' [`fabzB'] \\" _n
file write tab "\$c' + ab\$ & \textcolor{DeepNavy}{\bfseries `fsumA'} & \textcolor{DeepNavy}{\bfseries `fsumB'} \\" _n
file write tab "\$c\$: \$Y\$ on \$X\$ (total) & \textcolor{Maroon}{\bfseries `fcA'} (`fcAt') & \textcolor{Maroon}{\bfseries `fcB'} (`fcBt') \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_mediation_moderation_identity.tex

di "Done: tables/tab_mediation_moderation_identity.tex"


********************************************************************************
* MODULE 2: the sampling distribution of ab is skewed; Sobel versus percentile bootstrap intervals
********************************************************************************
********************************************************************************
* ACCTG 597A -- Week 6 (mediation and moderation module)
* The indirect effect ab is a product of two coefficients, so its sampling
* distribution is skewed in samples of the size common in practice. A
* normal-theory (Sobel) interval is symmetric around ab-hat and covers the
* true value less often than its nominal rate; a percentile bootstrap interval
* follows the shape of the sampling distribution.
*
* Part 1: Monte Carlo, 5,000 samples of N = 100 from
*   M = 0.25 X + e_M,  Y = 0.10 X + 0.25 M + e_Y,  all shocks standard normal,
* so ab = 0.0625. In every sample: the Sobel interval ab-hat +/- 1.96 s.e.
* and a percentile bootstrap interval from 1,000 resamples (computed in Mata
* for speed). Records the skewness of ab-hat and the coverage and rejection
* rate of both intervals. The figure draws the Monte Carlo distribution of
* ab-hat with the normal density implied by the mean Sobel s.e.
* Part 2: one further sample; the same bootstrap with Stata's bootstrap
* command, which is the recipe a researcher would use on real data.
* Checks: the Monte Carlo mean of ab-hat is within three Monte Carlo standard
* errors of ab; the bootstrap command's point estimate equals the sample ab.
* Requirements: base Stata. Output: figures/fig2_mediation_moderation_bootstrap.png
********************************************************************************

clear all
set seed 2017

local R  = 5000
local B  = 1000
local N  = 100
local a  = 0.25
local b  = 0.25
local cp = 0.10
local ab = `a'*`b'

* ---- Percentile bootstrap of ab, written in Mata so that 5,000 x 1,000 fits in minutes ----
mata:
real rowvector boot_ab(real colvector x, real colvector m, real colvector y,
                       real scalar B)
{
    real scalar    N, r, k
    real colvector idx, xb, mb, yb, ab, ba, bb
    real matrix    X1, X2
    N  = rows(x)
    ab = J(B, 1, .)
    for (r = 1; r <= B; r++) {
        idx = ceil(N :* runiform(N, 1))          // resample rows with replacement
        xb = x[idx]; mb = m[idx]; yb = y[idx]
        X1 = (J(N, 1, 1), xb)                     // M on X
        ba = invsym(cross(X1, X1)) * cross(X1, mb)
        X2 = (J(N, 1, 1), xb, mb)                 // Y on X and M
        bb = invsym(cross(X2, X2)) * cross(X2, yb)
        ab[r] = ba[2] * bb[3]
    }
    _sort(ab, 1)
    k = 0.025 * B                                 // 2.5th and 97.5th percentiles
    return(((ab[k] + ab[k + 1]) / 2, (ab[B - k] + ab[B - k + 1]) / 2))
}
end

* ---- Part 1: the sampling distribution of ab-hat and two intervals ----
tempname mc
postfile `mc' rep a_hat b_hat ab_hat sobel bl bu using mc_ab, replace
forvalues r = 1/`R' {
    quietly {
        clear
        set obs `N'
        gen x = rnormal()
        gen m = `a'*x + rnormal()
        gen y = `cp'*x + `b'*m + rnormal()
        reg m x
        local ah  = _b[x]
        local ase = _se[x]
        reg y x m
        local bh  = _b[m]
        local bse = _se[m]
        mata: st_matrix("ci", boot_ab(st_data(., "x"), st_data(., "m"), st_data(., "y"), `B'))
        post `mc' (`r') (`ah') (`bh') ((`ah')*(`bh')) (sqrt((`bh')^2*(`ase')^2 + (`ah')^2*(`bse')^2)) (ci[1,1]) (ci[1,2])
    }
}
postclose `mc'

use mc_ab, clear
count if missing(sobel) | missing(bl) | missing(bu)
di "Samples with a missing interval (should be 0): " r(N)
assert r(N) == 0
gen cov_sobel = (ab_hat - 1.96*sobel <= `ab') & (`ab' <= ab_hat + 1.96*sobel)
gen cov_boot  = (bl <= `ab') & (`ab' <= bu)
gen rej_sobel = abs(ab_hat/sobel) > 1.96
gen rej_boot  = (bl > 0) | (bu < 0)
sum ab_hat, detail
scalar mc_mean = r(mean)
scalar mc_sd   = r(sd)
scalar mc_skew = r(skewness)
sum sobel
scalar sobel_mean = r(mean)
sum cov_sobel
scalar S_cov = r(mean)
sum cov_boot
scalar B_cov = r(mean)
sum rej_sobel
scalar S_rej = r(mean)
sum rej_boot
scalar B_rej = r(mean)

di _n "Monte Carlo, `R' samples of N = `N': true ab = " %7.4f `ab'
di    "  mean ab-hat        = " %7.4f mc_mean "   (Monte Carlo s.e. of the mean " %7.4f mc_sd/sqrt(`R') ")"
di    "  s.d. of ab-hat     = " %7.4f mc_sd "   mean Sobel s.e. = " %7.4f sobel_mean
di    "  skewness of ab-hat = " %6.3f mc_skew "   (a normal distribution has skewness 0)"
di    "  95 percent interval coverage:  Sobel " %6.3f S_cov "   percentile bootstrap " %6.3f B_cov
di    "  rejection rate (interval excludes 0):  Sobel " %6.3f S_rej "   percentile bootstrap " %6.3f B_rej
di _n "Check: Monte Carlo mean within 3 s.e. of ab: " %8.5f abs(mc_mean - `ab') " < " %8.5f 3*mc_sd/sqrt(`R')
assert abs(mc_mean - `ab') < 3*mc_sd/sqrt(`R')

* ---- Figure: the Monte Carlo distribution of ab-hat against the normal ----
local fsk  : di %4.2f mc_skew
local fcs  : di %5.3f S_cov
local fcb  : di %5.3f B_cov
local fsd  : di %5.3f mc_sd
foreach s in fsk fcs fcb fsd {
    local `s' = strtrim("``s''")
}
local ytop = 1.20/(sobel_mean*sqrt(2*_pi))
capture mkdir figures
twoway (histogram ab_hat, width(0.008) color(navy%55) lcolor(white) lwidth(vthin))               ///
       (function normalden(x, `ab', sobel_mean), range(-0.08 0.22) lcolor(maroon) lwidth(thick)), ///
    xline(`ab', lcolor(gs8) lwidth(medthick))                                                     ///
    xline(0, lcolor(gs12) lpattern(dash))                                                         ///
    text(`=0.97*`ytop'' 0.152 "Monte Carlo distribution of ab-hat" "(`R' samples): skewness `fsk'", ///
         color(navy) size(medsmall) placement(e) justification(left))                             ///
    text(`=0.78*`ytop'' 0.152 "normal density with the mean" "Sobel s.e. (symmetric around ab)",  ///
         color(maroon) size(medsmall) placement(e) justification(left))                           ///
    text(`=0.57*`ytop'' 0.152 "95 percent interval coverage:" "Sobel `fcs', percentile bootstrap `fcb'", ///
         color(gs4) size(medsmall) placement(e) justification(left))                              ///
    xtitle("Indirect effect ab-hat") ytitle("Density")                                            ///
    xlabel(-0.05(0.05)0.20, format(%4.2f)) xscale(range(-0.08 0.30))                              ///
    ylabel(, angle(horizontal) nogrid) yscale(range(0 `=1.05*`ytop''))                            ///
    title("The sampling distribution of ab is skewed")                                            ///
    subtitle("5,000 simulated samples of N = 100 (seed 2017); solid vertical line at the true ab = 0.0625", size(small)) ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig2_mediation_moderation_bootstrap.png, replace width(2000)

* ---- Part 2: the recipe on one sample, with Stata's bootstrap command ----
clear
set obs `N'
gen x = rnormal()
gen m = `a'*x + rnormal()
gen y = `cp'*x + `b'*m + rnormal()
reg m x
scalar a1   = _b[x]
scalar a1se = _se[x]
reg y x m
scalar b1   = _b[m]
scalar b1se = _se[m]
scalar ab1  = a1*b1
scalar sob1 = sqrt(b1^2*a1se^2 + a1^2*b1se^2)

capture program drop indirect
program define indirect, rclass
    reg m x
    local ah = _b[x]
    reg y x m
    return scalar ab = `ah'*_b[m]
end
bootstrap ab = r(ab), reps(`B') saving(bs_ab, replace) nowarn nodots: indirect
estat bootstrap, percentile
scalar ab_bs = _b[ab]
di _n "Check: bootstrap point estimate minus sample ab = " %12.8f (ab_bs - ab1)
assert abs(ab_bs - ab1) < 1e-6
use bs_ab, clear
_pctile ab, p(2.5 97.5)
di _n "One sample, N = `N': ab = " %7.4f ab1 "   Sobel interval [" %7.4f ab1 - 1.96*sob1 ", " %7.4f ab1 + 1.96*sob1 "]" ///
      "   percentile bootstrap interval [" %7.4f r(r1) ", " %7.4f r(r2) "]"

capture erase mc_ab.dta
capture erase bs_ab.dta
di "Done: figures/fig2_mediation_moderation_bootstrap.png"


********************************************************************************
* MODULE 3: the conditional effect b1 + b3 W, centering, Johnson-Neyman boundaries, and a median split
********************************************************************************
********************************************************************************
* ACCTG 597A -- Week 6 (mediation and moderation module)
* In a moderation model Y = b1 X + b2 W + b3 XW + e, the effect of X on Y is
* the conditional effect theta(W) = b1 + b3 W: b1 is the effect at W = 0, and
* the Johnson-Neyman boundaries are the values of W where the conditional
* effect's t statistic equals the critical value.
*
* Pure simulation, N = 500, X and W standard normal:
*   Y = 0.15 X + 0.30 W + 0.25 X W + e,   e standard normal,
* so theta(W) = 0.15 + 0.25 W, zero at W = -0.6.
* Checks: mean-centered b1 equals b1 + b3 Wbar and b3 is unchanged (1e-6);
* at each Johnson-Neyman root, absolute theta divided by its s.e. equals the
* critical t (1e-6); lincom at W = 0 reproduces b1.
* Also compares a median split of W (two subgroup slopes) with the b3 test.
* Requirements: base Stata. Output: figures/fig3_mediation_moderation_conditional.png,
*   tables/tab_mediation_moderation_probe.tex
********************************************************************************

clear all
set seed 2017

set obs 500
gen x  = rnormal()
gen w  = rnormal()
gen xw = x*w
gen y  = 0.15*x + 0.30*w + 0.25*xw + rnormal()

* ---- The moderation regression ----
reg y x w xw
scalar b1  = _b[x]
scalar b2  = _b[w]
scalar b3  = _b[xw]
scalar b1t = _b[x]/_se[x]
scalar b3t = _b[xw]/_se[xw]
scalar tcrit = invttail(e(df_r), 0.025)
matrix V = e(V)
scalar v11 = V[1,1]
scalar v33 = V[3,3]
scalar v13 = V[1,3]

* ---- Pick-a-point: theta(W) at W = -1, 0, +1 (standard deviations of W) ----
foreach k in -1 0 1 {
    local kk = cond(`k' < 0, "m1", cond(`k' == 0, "0", "p1"))
    lincom x + (`k')*xw
    scalar th`kk'  = r(estimate)
    scalar th`kk't = r(estimate)/r(se)
}
di _n "Check: lincom at W = 0 minus b1 = " %12.8f (th0 - b1)
assert abs(th0 - b1) < 1e-6

* ---- Mean centering relabels b1 and leaves b3 unchanged ----
sum w
scalar wbar = r(mean)
gen wc  = w - wbar
gen xwc = x*wc
reg y x wc xwc
scalar b1c = _b[x]
scalar b3c = _b[xwc]
di "Check: centered b1 - (b1 + b3 Wbar) = " %12.8f (b1c - (b1 + b3*wbar))
di "Check: b3 - centered b3            = " %12.8f (b3 - b3c)
assert abs(b1c - (b1 + b3*wbar)) < 1e-6
assert abs(b3 - b3c) < 1e-6

* ---- Johnson-Neyman boundaries: (b1 + b3 W)^2 = t^2 Var(b1 + b3 W) ----
scalar qA = b3^2 - tcrit^2*v33
scalar qB = 2*(b1*b3 - tcrit^2*v13)
scalar qC = b1^2 - tcrit^2*v11
scalar disc = qB^2 - 4*qA*qC
assert disc > 0
scalar wjn1 = (-qB - sqrt(disc))/(2*qA)
scalar wjn2 = (-qB + sqrt(disc))/(2*qA)
scalar wlo = min(wjn1, wjn2)
scalar whi = max(wjn1, wjn2)
foreach r in lo hi {
    scalar th_`r' = b1 + b3*w`r'
    scalar se_`r' = sqrt(v11 + 2*w`r'*v13 + w`r'^2*v33)
    di "Check: at W = " %6.3f w`r' ", absolute t minus critical t = " %12.8f (abs(th_`r')/se_`r' - tcrit)
    assert abs(abs(th_`r')/se_`r' - tcrit) < 1e-6
}
di _n "b1 = " %6.3f b1 " (t = " %4.1f b1t "), b3 = " %6.3f b3 " (t = " %4.1f b3t "), critical t = " %6.4f tcrit
di    "theta(-1) = " %6.3f thm1 " (t = " %4.1f thm1t "); theta(0) = " %6.3f th0 " (t = " %4.1f th0t "); theta(+1) = " %6.3f thp1 " (t = " %4.1f thp1t ")"
di    "Johnson-Neyman boundaries: W = " %6.3f wlo " and W = " %6.3f whi
di    "  theta(W) differs from zero for W < " %6.3f wlo " (negative) and W > " %6.3f whi " (positive)"

* ---- A median split of W is not a test of b3 ----
sum w, detail
scalar wmed = r(p50)
reg y x if w <= wmed
scalar s_lo   = _b[x]
scalar s_lot  = _b[x]/_se[x]
scalar s_lose = _se[x]
reg y x if w > wmed
scalar s_hi   = _b[x]
scalar s_hit  = _b[x]/_se[x]
scalar s_hise = _se[x]
scalar dz = (s_hi - s_lo)/sqrt(s_hise^2 + s_lose^2)
di _n "Median split: slope of X below median W = " %6.3f s_lo " (t = " %4.1f s_lot ")"
di    "              slope of X above median W = " %6.3f s_hi " (t = " %4.1f s_hit ")"
di    "              difference z = " %4.1f dz "   versus b3 t = " %4.1f b3t

* ---- Table ----
foreach s in thm1 th0 thp1 wlo whi s_lo s_hi b3 {
    local f`s' : di %4.2f `s'
    local f`s' = strtrim("`f`s''")
}
foreach s in thm1t th0t thp1t s_lot s_hit dz b3t {
    local f`s' : di %3.1f `s'
    local f`s' = strtrim("`f`s''")
}
capture mkdir tables
file open tab using tables/tab_mediation_moderation_probe.tex, write replace
file write tab "\begin{tabular}{@{}lcc@{}}" _n
file write tab "\toprule" _n
file write tab " & Effect of \$X\$ & \$t\$ \\" _n
file write tab "\midrule" _n
file write tab "\$\theta(W)\$ at \$W = -1\$ & `fthm1' & `fthm1t' \\" _n
file write tab "\$\theta(W)\$ at \$W = 0\$ (\$= b_1\$) & `fth0' & `fth0t' \\" _n
file write tab "\$\theta(W)\$ at \$W = +1\$ & `fthp1' & `fthp1t' \\" _n
file write tab "Johnson-Neyman boundaries & \textcolor{Maroon}{\bfseries \$W = `fwlo'\$, \$`fwhi'\$} & \\" _n
file write tab "\midrule" _n
file write tab "Median split: \$W\$ below median & `fs_lo' & `fs_lot' \\" _n
file write tab "Median split: \$W\$ above median & `fs_hi' & `fs_hit' \\" _n
file write tab "Difference of subgroup slopes & & `fdz' \\" _n
file write tab "\$b_3\$ (interaction) & \textcolor{DeepNavy}{\bfseries `fb3'} & \textcolor{DeepNavy}{\bfseries `fb3t'} \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_mediation_moderation_probe.tex

* ---- Figure: the conditional effect with its band and the JN boundaries ----
local fb1 : di %4.2f b1
local fb3 : di %4.2f b3
local fb1 = strtrim("`fb1'")
local fb3 = strtrim("`fb3'")
clear
set obs 101
gen s  = -2.5 + 5*(_n - 1)/100
gen me = b1 + b3*s
gen se = sqrt(v11 + 2*s*v13 + s^2*v33)
gen lo = me - tcrit*se
gen hi = me + tcrit*se
gen mid = (s >= wlo & s <= whi)

capture mkdir figures
twoway (rarea lo hi s if !mid, color(navy%15) lwidth(none))                           ///
       (rarea lo hi s if mid, color(gs14%60) lwidth(none))                             ///
       (line me s, lcolor(navy) lwidth(thick)),                                        ///
    yline(0, lcolor(gs8) lpattern(dash))                                               ///
    xline(`=wlo' `=whi', lcolor(maroon) lpattern(shortdash) lwidth(medthick))          ///
    text(-0.45 0.3 "{&theta}(W) = b{sub:1} + b{sub:3} W = `fb1' + `fb3' W", color(navy) size(medium) placement(e)) ///
    text(0.62 `=(wlo+whi)/2' "{&theta}(W) not" "different" "from zero", color(gs5) size(small)) ///
    text(-0.72 `=wlo-0.08' "W = `fwlo'", color(maroon) size(medsmall) placement(w))     ///
    text(-0.72 `=whi+0.08' "W = `fwhi'", color(maroon) size(medsmall) placement(e))     ///
    text(-0.98 `=whi+0.08' "Johnson-Neyman boundaries", color(maroon) size(medsmall) placement(e)) ///
    xtitle("Moderator W (standard deviations)") ytitle("Conditional effect of X on Y")  ///
    xlabel(-2(1)2) ylabel(-1(0.5)1, angle(horizontal) format(%3.1f))                    ///
    yscale(range(-1.1 1.1))                                                            ///
    title("The effect of X on Y is a function of W")                                   ///
    subtitle("Simulated data (seed 2017), N = 500; band is a 95 percent confidence interval") ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig3_mediation_moderation_conditional.png, replace width(2000)

di "Done: figures/fig3_mediation_moderation_conditional.png"


********************************************************************************
* MODULE 4: moderated mediation, the conditional indirect effect, and the index a3 b
********************************************************************************
********************************************************************************
* ACCTG 597A -- Week 6 (mediation and moderation module)
* Moderated mediation: when W moderates the a path, the indirect effect of X
* on Y through M is a function of W, omega(W) = (a1 + a3 W) b, and the index
* of moderated mediation a3 b is its slope.
*
* Pure simulation, N = 500, X and W standard normal:
*   M = 0.5 X + 0.2 W + 0.3 X W + e_M
*   Y = 0.2 X + 0.6 M + e_Y                 (W and XW enter the estimated Y
*                                            equation with true coefficients 0)
* so omega(W) = 0.30 + 0.18 W and the index is 0.18.
* Checks (1e-6): the reduced-form regression of Y on X, W, XW returns
* c1' + b a1 on X and c3' + b a3 on XW; omega(1) - omega(0) equals a3 b.
* Inference: 1,000 bootstrap resamples of (a1, a3, b); percentile band for
* omega(W) on a grid and a percentile interval for the index.
* Requirements: base Stata. Output: figures/fig4_mediation_moderation_modmed.png
********************************************************************************

clear all
set seed 2017

set obs 500
gen x  = rnormal()
gen w  = rnormal()
gen xw = x*w
gen m  = 0.5*x + 0.2*w + 0.3*xw + rnormal()
gen y  = 0.2*x + 0.6*m + rnormal()

* ---- The two path equations (PROCESS model 8 of Hayes 2013: W may moderate both a and c') ----
reg m x w xw
scalar a1 = _b[x]
scalar a2 = _b[w]
scalar a3 = _b[xw]
scalar a3t = _b[xw]/_se[xw]
reg y x w xw m
scalar c1 = _b[x]
scalar c2 = _b[w]
scalar c3 = _b[xw]
scalar b  = _b[m]
scalar bt = _b[m]/_se[m]
scalar om0   = a1*b
scalar index = a3*b

* ---- Identity checks against the reduced form ----
reg y x w xw
di _n "Identity checks (should be zero to machine precision):"
di "  reduced-form X  minus (c1' + b a1) = " %12.8f (_b[x] - (c1 + b*a1))
di "  reduced-form XW minus (c3' + b a3) = " %12.8f (_b[xw] - (c3 + b*a3))
di "  reduced-form W  minus (c2' + b a2) = " %12.8f (_b[w] - (c2 + b*a2))
assert abs(_b[x]  - (c1 + b*a1)) < 1e-6
assert abs(_b[xw] - (c3 + b*a3)) < 1e-6
assert abs(_b[w]  - (c2 + b*a2)) < 1e-6
di "  omega(1) - omega(0) minus a3 b     = " %12.8f (((a1 + a3)*b - a1*b) - index)
assert abs(((a1 + a3)*b - a1*b) - index) < 1e-6

di _n "a1 = " %6.3f a1 ", a3 = " %6.3f a3 " (t = " %4.1f a3t "), b = " %6.3f b " (t = " %4.1f bt ")"
di    "conditional indirect effect at W = 0: a1 b = " %6.3f om0
di    "index of moderated mediation a3 b     = " %6.3f index
di    "conditional indirect effect at W = -1, 0, +1: " %6.3f (a1 - a3)*b ", " %6.3f om0 ", " %6.3f (a1 + a3)*b

* ---- Bootstrap the three path coefficients ----
capture program drop modmed
program define modmed, rclass
    reg m x w xw
    return scalar a1 = _b[x]
    return scalar a3 = _b[xw]
    reg y x w xw m
    return scalar b  = _b[m]
end
bootstrap a1 = r(a1) a3 = r(a3) b = r(b), reps(1000) saving(bs_mm, replace) nowarn nodots: modmed
di _n "Check: bootstrap point estimates equal the sample estimates: " ///
      %12.8f (_b[a1] - a1) " " %12.8f (_b[a3] - a3) " " %12.8f (_b[b] - b)
assert abs(_b[a1] - a1) < 1e-6 & abs(_b[a3] - a3) < 1e-6 & abs(_b[b] - b) < 1e-6

use bs_mm, clear
rename (a1 a3 b) (a1_bs a3_bs b_bs)   // avoid shadowing the scalars a1, a3, b
gen idx_bs = a3_bs*b_bs
_pctile idx_bs, p(2.5 97.5)
scalar ixl = r(r1)
scalar ixu = r(r2)
di "index of moderated mediation: " %6.3f index "   percentile 95 percent interval [" %6.3f ixl ", " %6.3f ixu "]"

* percentile band for omega(W) on a grid
tempname grid
postfile `grid' s om lo hi using om_grid, replace
forvalues k = 0/40 {
    local s = -2 + 4*`k'/40
    quietly {
        gen om_k = (a1_bs + a3_bs*`s')*b_bs
        _pctile om_k, p(2.5 97.5)
        post `grid' (`s') ((a1 + a3*`s')*b) (r(r1)) (r(r2))
        drop om_k
    }
}
postclose `grid'
use om_grid, clear
gen truth = 0.30 + 0.18*s

* ---- Figure ----
local fidx : di %4.2f index
local fixl : di %4.2f ixl
local fixu : di %4.2f ixu
local fom0 : di %4.2f om0
foreach s in fidx fixl fixu fom0 {
    local `s' = strtrim("``s''")
}
capture mkdir figures
twoway (rarea lo hi s, color(navy%15) lwidth(none))                                      ///
       (line om s, lcolor(navy) lwidth(thick))                                          ///
       (line truth s, lcolor(maroon) lpattern(dash) lwidth(medthick)),                  ///
    yline(0, lcolor(gs8) lpattern(dash))                                                ///
    text(0.86 -1.95 "{&omega}(W) = (a{sub:1} + a{sub:3} W) b: estimate with 95 percent band", color(navy) size(medium) placement(e)) ///
    text(0.74 -1.95 "index of moderated mediation a{sub:3} b = `fidx' [`fixl', `fixu']", color(navy) size(medium) placement(e)) ///
    text(-0.22 1.95 "true {&omega}(W) = 0.30 + 0.18 W", color(maroon) size(medium) placement(w)) ///
    xtitle("Moderator W (standard deviations)") ytitle("Indirect effect of X through M")   ///
    xlabel(-2(1)2) ylabel(-0.4(0.2)1, angle(horizontal) format(%3.1f))                   ///
    yscale(range(-0.45 0.95))                                                           ///
    title("The indirect effect through M rises with W")                                 ///
    subtitle("Simulated data (seed 2017), N = 500; band is a 95 percent percentile bootstrap interval (1,000 resamples)", size(small)) ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig4_mediation_moderation_modmed.png, replace width(2000)

capture erase bs_mm.dta
capture erase om_grid.dta
di "Done: figures/fig4_mediation_moderation_modmed.png"


********************************************************************************
* MODULE 5: the identity fixes the sum, not the split: a shared unobservable and measurement error in M
********************************************************************************
********************************************************************************
* ACCTG 597A -- Week 6 (mediation and moderation module)
* The decomposition c = c' + ab holds in every sample by construction; what
* moves is the split. An unobserved variable in both M and Y, or measurement
* error in M, reallocates the total between the direct and indirect paths
* while the total effect of an exogenous X stays unbiased.
*
* Monte Carlo, 500 samples of N = 2,000; X, U, all shocks standard normal.
*   Design A (clean):        M = 0.5 X + e_M,           Y = 0.2 X + 0.8 M + e_Y
*   Design B (U in M and Y): M = 0.5 X + 0.6 U + e_M,   Y = 0.2 X + 0.8 M + 0.6 U + e_Y
*   Design C (noisy M):      observed M* = M_A + u,      Y = Y_A
* True direct c' = 0.2, indirect ab = 0.4, total c = 0.6 in every design.
* Population OLS values (derived in mediation_moderation_derivations.py):
*   B: b -> 0.8 + 0.36/1.36 = 1.0647, ab -> 0.5324, c' -> 0.0676
*   C: b -> 0.8 x 0.5 = 0.40,         ab -> 0.2000, c' -> 0.4000
* Checks: c - (c' + ab) = 0 (1e-6) in every sample and design; Monte Carlo
* means within 0.01 of the population values.
* Requirements: base Stata. Output: figures/fig5_mediation_moderation_paths.png
********************************************************************************

clear all
set seed 2012

local R = 500
local N = 2000

tempname pf
postfile `pf' rep dA iA cA dB iB cB dC iC cC viol using paths_mc, replace
forvalues r = 1/`R' {
    quietly {
        clear
        set obs `N'
        gen x  = rnormal()
        gen u  = rnormal()
        gen em = rnormal()
        gen ey = rnormal()
        gen mn = rnormal()
        gen mA = 0.5*x + em
        gen yA = 0.2*x + 0.8*mA + ey
        gen mB = 0.5*x + 0.6*u + em
        gen yB = 0.2*x + 0.8*mB + 0.6*u + ey
        gen mC = mA + mn
        gen yC = yA
        local viol = 0
        foreach d in A B C {
            reg m`d' x
            local a  = _b[x]
            reg y`d' x m`d'
            local cp = _b[x]
            local b  = _b[m`d']
            reg y`d' x
            local c  = _b[x]
            local d`d' = `cp'
            local i`d' = `a'*`b'
            local c`d' = `c'
            if abs(`c' - (`cp' + `a'*`b')) > 1e-6 local viol = `viol' + 1
        }
        post `pf' (`r') (`dA') (`iA') (`cA') (`dB') (`iB') (`cB') (`dC') (`iC') (`cC') (`viol')
    }
}
postclose `pf'

use paths_mc, clear
sum viol
di _n "Identity check: samples in which c - (c' + ab) exceeded 1e-6 (any design): " r(max)
assert r(max) == 0
collapse (mean) dA iA cA dB iB cB dC iC cC
foreach v in dA iA cA dB iB cB dC iC cC {
    scalar m_`v' = `v'[1]
}

* population values
scalar biasB = 0.6*0.6/(0.36 + 1)
scalar pop_iB = 0.5*(0.8 + biasB)
scalar pop_dB = 0.2 - 0.5*biasB
scalar lam    = 1/(1 + 1)
scalar pop_iC = 0.4*lam
scalar pop_dC = 0.2 + 0.4*(1 - lam)

di _n "Monte Carlo means (`R' samples of N = `N'):"
di    "  Design A: direct " %6.3f m_dA "  indirect " %6.3f m_iA "  total " %6.3f m_cA "   (true 0.200, 0.400, 0.600)"
di    "  Design B: direct " %6.3f m_dB "  indirect " %6.3f m_iB "  total " %6.3f m_cB "   (population " %6.3f pop_dB ", " %6.3f pop_iB ", 0.600)"
di    "  Design C: direct " %6.3f m_dC "  indirect " %6.3f m_iC "  total " %6.3f m_cC "   (population " %6.3f pop_dC ", " %6.3f pop_iC ", 0.600)"
di _n "Checks (Monte Carlo mean minus population value, tolerance 0.01):"
di    "  A direct " %8.4f (m_dA - 0.2) "  A indirect " %8.4f (m_iA - 0.4)
di    "  B direct " %8.4f (m_dB - pop_dB) "  B indirect " %8.4f (m_iB - pop_iB)
di    "  C direct " %8.4f (m_dC - pop_dC) "  C indirect " %8.4f (m_iC - pop_iC)
di    "  totals   " %8.4f (m_cA - 0.6) " " %8.4f (m_cB - 0.6) " " %8.4f (m_cC - 0.6)
assert abs(m_dA - 0.2) < 0.01 & abs(m_iA - 0.4) < 0.01
assert abs(m_dB - pop_dB) < 0.01 & abs(m_iB - pop_iB) < 0.01
assert abs(m_dC - pop_dC) < 0.01 & abs(m_iC - pop_iC) < 0.01
assert abs(m_cA - 0.6) < 0.01 & abs(m_cB - 0.6) < 0.01 & abs(m_cC - 0.6) < 0.01

* ---- Figure: stacked bars, direct (maroon) below indirect (navy) ----
clear
set obs 3
gen pos = _n
gen direct = .
gen indirect = .
replace direct = m_dA in 1
replace direct = m_dB in 2
replace direct = m_dC in 3
replace indirect = m_iA in 1
replace indirect = m_iB in 2
replace indirect = m_iC in 3
gen total = direct + indirect
gen ymid_d = direct/2
gen ymid_i = direct + indirect/2
gen lab_d = string(direct, "%4.2f")
gen lab_i = string(indirect, "%4.2f")
gen share = string(100*indirect/total, "%2.0f") + " percent mediated"

capture mkdir figures
twoway (bar direct pos, barwidth(0.55) color(maroon))                                    ///
       (rbar direct total pos, barwidth(0.55) color(navy))                               ///
       (scatter ymid_d pos, msymbol(none) mlabel(lab_d) mlabposition(0) mlabcolor(white) mlabsize(medium)) ///
       (scatter ymid_i pos, msymbol(none) mlabel(lab_i) mlabposition(0) mlabcolor(white) mlabsize(medium)) ///
       (scatter total pos, msymbol(none) mlabel(share) mlabposition(12) mlabcolor(gs4) mlabsize(medsmall)), ///
    yline(0.6, lcolor(gs8) lpattern(dash))                                              ///
    yline(0.2, lcolor(maroon) lpattern(shortdash))                                      ///
    text(0.572 3.32 "total c = 0.60", color(gs4) size(medsmall) placement(e))                 ///
    text(0.50 3.32 "indirect ab", color(navy) size(medium) placement(e))                     ///
    text(0.165 3.32 "true direct c' = 0.20", color(maroon) size(medsmall) placement(e))      ///
    text(0.09 3.32 "direct c'", color(maroon) size(medium) placement(e))                     ///
    xlabel(1 `""A: clean" "design""' 2 `""B: U enters" "both M and Y""' 3 `""C: M measured" "with error""', noticks) ///
    xscale(range(0.5 4.0))                                                              ///
    ylabel(0(0.2)0.8, angle(horizontal) format(%3.1f)) yscale(range(0 0.8))             ///
    xtitle("") ytitle("Estimated effect of X on Y")                                     ///
    title("The split between paths moves; the total does not")                          ///
    subtitle("Monte Carlo means, 500 samples of N = 2,000 (seed 2012); true c' = 0.20, ab = 0.40, c = 0.60") ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig5_mediation_moderation_paths.png, replace width(2000)

capture erase paths_mc.dta
di "Done: figures/fig5_mediation_moderation_paths.png"


log close
