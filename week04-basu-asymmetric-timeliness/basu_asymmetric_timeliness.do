********************************************************************************
* ACCTG 597A -- Week 4: Reversal and Truncation
* Basu (1997), "The Conservatism Principle and the Asymmetric Timeliness of
* Earnings," Journal of Accounting and Economics 24(1): 3-37, read with
* Dietrich, Muller, and Riedl (2007), "Asymmetric Timeliness Tests of
* Accounting Conservatism," Review of Accounting Studies 12(1): 95-124 (DMR).
*
* Self-contained simulation, base Stata only (version 17 or later). Run from
* this folder:  do basu_asymmetric_timeliness.do
* Figures are written to figures/ and booktabs table fragments to tables/.
* Runtime about three minutes (Module 5 is a 7,000-sample Monte Carlo).
*
* Modules. In each, the true parameter is known by construction (returns
* R = X + eta with beta = lambda = 1 and cov(X, eta) = 0, so earnings are not
* conservative and the true asymmetric timeliness coefficient is zero), so
* the property under study is visible directly. Seed 1997 throughout.
*   1. Identity. N = 50,000 with DMR Table 1 moments; each sign-partitioned
*      reverse slope equals the full-sample forward slope + SVR + ST exactly;
*      Basu's beta1 equals the difference of the two subsample slopes.
*   2. Truncation geometry. Unit normals; splitting on R = X + eta makes
*      cov(X, eta) = -1/pi in both halves although it is zero overall.
*   3. DMR Fig. 3 examples. The three printed twelve-observation panels
*      reproduced to 0.001, plus the mirror of Panel C (sign flips exactly).
*   4. DMR Table 2. The placebo-return columns recreated with the
*      full + SVR + ST arithmetic asserted.
*   5. Monte Carlo. Eight designs, 1,000 samples of N = 2,000: normal data
*      give a coefficient centered at zero whether the split is at or away
*      from the mean; the six non-normal designs range from about -0.26 to
*      +0.82 with rejection rates of 66 to 100 percent.
********************************************************************************

version 17
clear all
set more off
capture log close
log using basu_asymmetric_timeliness.log, replace text
capture mkdir figures
capture mkdir tables

********************************************************************************
* MODULE 1: The exact sample decomposition of the two Basu slopes
********************************************************************************
* Week 4, Basu (1997): the exact sample decomposition of the two Basu slopes
*
* Demonstrates that within each sign-partitioned subsample the reverse slope
* of earnings on returns equals the full-sample forward slope plus a sample-
* variance-ratio (SVR) term plus a sample-truncation (ST) term, exactly, and
* that Basu's interaction coefficient is the difference of the two subsample
* slopes (Dietrich, Muller, and Riedl 2007, eqs. 1.7a-b and Table 2 Panels D-E).
*
* DGP (DMR eq. 1.1 with beta = lambda = 1, no conservatism): R = X + eta,
* cov(X, eta) = 0, N = 50,000, moments calibrated to DMR Table 1:
* sd(X) = 0.17 (var 0.029), sd(eta) = 0.497 (var 0.247, so var(R) = 0.276).
*   Design 1: X and eta symmetric normal, split at the mean of R (zero).
*   Design 2: X normal, eta right-skewed (standardized chi-square(2), mean 0,
*             sd 0.497, skewness 2), split at zero; returns then have the long
*             right tail of DMR Table 1.
* Identities checked (zero to machine precision):
*   full sample: delta = b s2_X/s2_R and R2 = delta b (DMR eq. 1.4; Basu fn. 7)
*   each half j: delta_j = b + SVR_j + ST_j, SVR_j = b s2_X,j/s2_R,j - b,
*                ST_j = s_Xtau,j/s2_R,j, tau the full-sample forward residual
*   each half j: delta_j^2 = R2_j s2_X,j/s2_R,j (Basu fn. 7)
*   stacked regression: beta1 = delta_1 - delta_0, beta0 = delta_0 (DMR fn. 14)
* Requirements: base Stata. Output: tables/tab_basu_identity.tex

clear all
set seed 1997

local N = 50000

forvalues d = 1/2 {
    clear
    set obs `N'
    gen X = 0.17*rnormal()
    if `d' == 1 gen eta = 0.497*rnormal()
    if `d' == 2 gen eta = 0.497*(rchi2(2) - 2)/2
    gen R  = X + eta
    gen DR = (R < 0)
    gen RDR = R*DR

    * full-sample forward and reverse regressions
    reg R X
    scalar b`d' = _b[X]
    scalar r2f`d' = e(r2)
    predict tau, resid
    reg X R
    scalar dfull`d' = _b[R]
    scalar r2r`d' = e(r2)
    qui corr X R eta tau, cov
    matrix C = r(C)
    scalar vXf`d' = C[1,1]
    scalar vRf`d' = C[2,2]
    scalar cXetaf`d' = C[1,3]

    * stacked Basu regression
    reg X DR R RDR
    scalar bas0`d' = _b[R]
    scalar bas1`d' = _b[RDR]
    scalar bas1t`d' = _b[RDR]/_se[RDR]

    * subsamples: j = 0 good news (R >= 0), j = 1 bad news (R < 0)
    forvalues j = 0/1 {
        if `j' == 0 local cond "R >= 0"
        else local cond "R < 0"
        reg X R if `cond'
        scalar del`d'`j' = _b[R]
        scalar r2`d'`j'  = e(r2)
        qui corr X R eta tau if `cond', cov
        matrix C = r(C)
        scalar vX`d'`j'    = C[1,1]
        scalar vR`d'`j'    = C[2,2]
        scalar cXeta`d'`j' = C[1,3]
        scalar cXtau`d'`j' = C[1,4]
        scalar svr`d'`j'   = b`d'*vX`d'`j'/vR`d'`j' - b`d'
        scalar st`d'`j'    = cXtau`d'`j'/vR`d'`j'
        scalar sum`d'`j'   = b`d' + svr`d'`j' + st`d'`j'
    }

    di _n "Design `d' identity checks (should be zero to machine precision):"
    di "  full: delta minus b*s2X/s2R          : " %10.6f (dfull`d' - b`d'*vXf`d'/vRf`d')
    di "  full: R2(reverse) minus R2(forward)  : " %10.6f (r2r`d' - r2f`d')
    di "  full: R2 minus delta*b               : " %10.6f (r2r`d' - dfull`d'*b`d')
    di "  full: cov(X, tau) (OLS orthogonality): " %10.6f C[1,4]
    forvalues j = 0/1 {
        di "  half `j': delta_j minus (b + SVR_j + ST_j): " %10.6f (del`d'`j' - sum`d'`j')
        di "  half `j': delta_j^2 minus R2_j*s2X_j/s2R_j : " %10.6f (del`d'`j'^2 - r2`d'`j'*vX`d'`j'/vR`d'`j')
        assert abs(del`d'`j' - sum`d'`j') < 1e-6
        assert abs(del`d'`j'^2 - r2`d'`j'*vX`d'`j'/vR`d'`j') < 1e-6
    }
    di "  stacked: beta1 minus (delta_1 - delta_0): " %10.6f (bas1`d' - (del`d'1 - del`d'0))
    di "  stacked: beta0 minus delta_0            : " %10.6f (bas0`d' - del`d'0)
    assert abs(dfull`d' - b`d'*vXf`d'/vRf`d') < 1e-6
    assert abs(r2r`d' - r2f`d') < 1e-6
    assert abs(r2r`d' - dfull`d'*b`d') < 1e-6
    assert abs(bas1`d' - (del`d'1 - del`d'0)) < 1e-6
    assert abs(bas0`d' - del`d'0) < 1e-6

    di _n "Design `d' values:"
    di "  forward slope b = " %6.3f b`d' "   full-sample reverse slope delta = " %6.3f dfull`d' "   s2X/s2R = " %6.3f vXf`d'/vRf`d'
    forvalues j = 0/1 {
        di "  half `j': delta_j = " %7.3f del`d'`j' "  SVR_j = " %7.3f svr`d'`j' "  ST_j = " %7.3f st`d'`j' ///
           "  cov(X,eta|j) = " %8.4f cXeta`d'`j' "  R2_j = " %6.3f r2`d'`j' "  s2X_j/s2R_j = " %6.3f vX`d'`j'/vR`d'`j'
    }
    di "  Basu beta1 = " %6.3f bas1`d' " (t = " %5.1f bas1t`d' ")"
}

* ---- numbers that travel to the slide: formats that do not pad ----
foreach s in b1 b2 dfull1 dfull2 bas11 bas12 {
    local f`s' : di %6.3f `s'
    local f`s' = strtrim("`f`s''")
}
forvalues d = 1/2 {
    forvalues j = 0/1 {
        foreach s in del svr st sum r2 {
            local f`s'`d'`j' : di %6.3f `s'`d'`j'
            local f`s'`d'`j' = strtrim("`f`s'`d'`j''")
        }
        local fcx`d'`j' : di %6.3f cXeta`d'`j'
        local fcx`d'`j' = strtrim("`fcx`d'`j''")
    }
    local fbt`d' : di %4.1f bas1t`d'
    local fbt`d' = strtrim("`fbt`d''")
}

* ---- booktabs fragment: escape math dollars as \$ (file write expands $) ----
file open tab using tables/tab_basu_identity.tex, write replace
file write tab "\begin{tabular}{@{}lcccc@{}}" _n
file write tab "\toprule" _n
file write tab " & \multicolumn{2}{c}{Symmetric normal} & \multicolumn{2}{c}{Right-skewed \$\eta\$} \\" _n
file write tab "\cmidrule(lr){2-3} \cmidrule(lr){4-5}" _n
file write tab " & Good, \$R \ge 0\$ & Bad, \$R < 0\$ & Good, \$R \ge 0\$ & Bad, \$R < 0\$ \\" _n
file write tab "\midrule" _n
file write tab "Forward slope \$b\$ & `fb1' & `fb1' & `fb2' & `fb2' \\" _n
file write tab "SVR\$_j\$ & `fsvr10' & `fsvr11' & `fsvr20' & `fsvr21' \\" _n
file write tab "ST\$_j\$ & `fst10' & `fst11' & `fst20' & `fst21' \\" _n
file write tab "\$b + \mathrm{SVR}_j + \mathrm{ST}_j\$ & `fsum10' & `fsum11' & `fsum20' & `fsum21' \\" _n
file write tab "Basu slope \$\hat\delta_j\$ & \textcolor{DeepNavy}{\bfseries `fdel10'} & \textcolor{DeepNavy}{\bfseries `fdel11'} & \textcolor{Maroon}{\bfseries `fdel20'} & \textcolor{Maroon}{\bfseries `fdel21'} \\" _n
file write tab "\midrule" _n
file write tab "\$\mathrm{cov}(X, \eta \mid j)\$ & `fcx10' & `fcx11' & `fcx20' & `fcx21' \\" _n
file write tab "\$R^2_j\$ & `fr210' & `fr211' & `fr220' & `fr221' \\" _n
file write tab "\midrule" _n
file write tab "\$\hat\beta_1\$ (\$t\$) & \multicolumn{2}{c}{\textcolor{DeepNavy}{\bfseries `fbas11' (`fbt1')}} & \multicolumn{2}{c}{\textcolor{Maroon}{\bfseries `fbas12' (`fbt2')}} \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_basu_identity.tex
di "Done: tables/tab_basu_identity.tex"


********************************************************************************
* MODULE 2: Sorting on the sign of R = X + eta tilts each half of the cloud
********************************************************************************
* Week 4, Basu (1997): sorting on the sign of R = X + eta tilts each half
*
* Demonstrates the sample-truncation mechanism of Dietrich, Muller, and Riedl
* (2007, Fig. 2, p. 103): when the sample is split on the sign of returns,
* which contain the non-earnings component eta, the covariance between X and
* eta is negative in both halves although it is zero in the full sample.
*
* DGP (stylized for the picture): X and eta independent standard normal,
* N = 4,000, R = X + eta, no conservatism. The split R < 0 is the half-plane
* eta < -X. Under unit normals the within-half covariance has the closed form
* cov(X, eta | R < 0) = cov(X, eta | R >= 0) = -1/pi = -0.318 (from
* R = X + eta and W = X - eta independent), which the script compares with the
* sample value at a tolerance of 0.04 (about two Monte Carlo standard errors).
* Requirements: base Stata. Output: figures/fig1_basu_truncation.png

clear all
set seed 1997

set obs 4000
gen X   = rnormal()
gen eta = rnormal()
gen R   = X + eta
gen bad = (R < 0)

* full-sample covariance and slope of eta on X (both near zero)
qui corr X eta, cov
matrix C = r(C)
scalar cfull = C[1,2]
reg eta X
scalar sfull = _b[X]

* within-half covariances, correlations, and slopes of eta on X
forvalues j = 0/1 {
    qui corr X eta if bad == `j', cov
    matrix C = r(C)
    scalar c`j' = C[1,2]
    scalar v`j' = C[1,1]
    qui corr X eta if bad == `j'
    scalar rho`j' = r(rho)
    reg eta X if bad == `j'
    scalar s`j' = _b[X]
    scalar a`j' = _b[_cons]
}

di _n "Truncation checks:"
di "  full-sample cov(X, eta)              : " %8.4f cfull
di "  cov(X, eta | R >= 0)                 : " %8.4f c0 "   analytical -1/pi = " %8.4f -1/_pi
di "  cov(X, eta | R < 0)                  : " %8.4f c1 "   analytical -1/pi = " %8.4f -1/_pi
di "  slope of eta on X, good / bad halves : " %8.4f s0 " / " %8.4f s1 "   analytical -1/(pi-1) = " %8.4f -1/(_pi-1)
di "  corr(X, eta | R >= 0) / (R < 0)      : " %8.4f rho0 " / " %8.4f rho1
assert abs(cfull) < 0.06
assert abs(c0 - (-1/_pi)) < 0.04
assert abs(c1 - (-1/_pi)) < 0.04

* binned means of eta by X within each half, for the picture
gen bin = .
forvalues j = 0/1 {
    xtile b`j' = X if bad == `j', nq(12)
    replace bin = b`j' if bad == `j'
}
bysort bad bin: egen mX = mean(X)
bysort bad bin: egen mE = mean(eta)
egen tag = tag(bad bin)

local fc0 : di %5.2f c0
local fc0 = strtrim("`fc0'")
local fc1 : di %5.2f c1
local fc1 = strtrim("`fc1'")
local fs0 : di %5.2f s0
local fs0 = strtrim("`fs0'")
local fs1 : di %5.2f s1
local fs1 = strtrim("`fs1'")
local fcf : di %5.2f cfull
local fcf = strtrim("`fcf'")

twoway (scatter eta X if bad == 0, mcolor(navy%18) msymbol(oh) msize(vsmall))          ///
       (scatter eta X if bad == 1, mcolor(maroon%18) msymbol(oh) msize(vsmall))        ///
       (function y = -x, range(-3.6 3.6) lcolor(gs8) lwidth(medium))                    ///
       (function y = a0 + s0*x, range(-2.6 3.4) lcolor(navy) lwidth(thick))             ///
       (function y = a1 + s1*x, range(-3.4 2.6) lcolor(maroon) lwidth(thick))           ///
       (scatter mE mX if tag & bad == 0, mcolor(navy) msymbol(circle) msize(medium))    ///
       (scatter mE mX if tag & bad == 1, mcolor(maroon) msymbol(circle) msize(medium)), ///
    yline(0, lcolor(gs10) lpattern(dash))                                               ///
    text(3.2 1.4 "good news, R {&ge} 0" "cov(X, {&eta}) = `fc0', slope `fs0'",           ///
         color(navy) size(medsmall) placement(e) justification(left))                   ///
    text(-3.2 -3.5 "bad news, R < 0" "cov(X, {&eta}) = `fc1', slope `fs1'",              ///
         color(maroon) size(medsmall) placement(e) justification(left))                 ///
    text(1.6 -3.5 "boundary {&eta} = -X",                                                 ///
         color(gs6) size(medsmall) placement(e) justification(left))                    ///
    text(3.55 -1.7 "full sample: cov(X, {&eta}) = `fcf'",                                 ///
         color(gs6) size(medsmall) placement(e) justification(left))                    ///
    xtitle("Earnings X") ytitle("Non-earnings information {&eta}")                     ///
    xlabel(-3(1)3) ylabel(-3(1)3, angle(horizontal) format(%2.0f))                      ///
    title("A split on the sign of R tilts both halves of the cloud")                   ///
    subtitle("Simulated data (seed 1997), N = 4,000; X, {&eta} independent standard normal; R = X + {&eta}") ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig1_basu_truncation.png, replace width(2000)
di "Done: figures/fig1_basu_truncation.png"


********************************************************************************
* MODULE 3: DMR's twelve-observation examples, recreated and mirrored
********************************************************************************
* Week 4, Basu (1997): the twelve-observation examples of Dietrich, Muller,
* and Riedl (2007, Fig. 3, pp. 105-111), recreated and extended
*
* Reproduces, from the printed data, the full-sample and subsample slopes of
* the three numerical examples (Panel A: X symmetric and eta homoskedastic;
* Panel B: eta variance rising in X; Panel C: X left-skewed) and adds the
* mirror image of Panel C (X and eta negated), in which the two subsamples
* swap and the Basu coefficient changes sign exactly.
*
* Data: R = X + eta with beta = lambda = 1 and cov(X, eta) = 0 in each panel.
* Checks: the recreated slopes match DMR's printed 0.429, 0.813, 0.813
* (Panel A), 0.038 and 1.917 (Panel B), 0.281 and 0.602 (Panel C) to 0.001;
* Panel A moments 13.636, 31.818, 13.367, 10.7, 1.249, 0.651 to 0.001; the
* mirror identities delta_1(mirror) = delta_0(C) and delta_0(mirror) =
* delta_1(C) hold to machine precision. The ST ratio is computed with n - 1
* for both covariance and variance (DMR fn. 17 print -0.363 with n for the
* covariance; with consistent scaling it is -0.436 and the identity is exact).
* Requirements: base Stata. Outputs: tables/tab_basu_dmr_a.tex,
* tables/tab_basu_dmr_bc.tex

clear all

capture program drop basupair
program define basupair, rclass
    * computes full-sample forward slope b, reverse slope delta, and the two
    * subsample reverse slopes delta0 (R >= 0) and delta1 (R < 0) from X, eta
    syntax varlist(min=2 max=2)
    tokenize `varlist'
    tempvar R tau
    gen `R' = `1' + `2'
    qui reg `R' `1'
    return scalar b = _b[`1']
    qui predict `tau', resid
    qui reg `1' `R'
    return scalar delta = _b[`R']
    qui corr `1' `R' `2' `tau', cov
    matrix C = r(C)
    return scalar vXf = C[1,1]
    return scalar vRf = C[2,2]
    forvalues j = 0/1 {
        if `j' == 0 local cond "`R' >= 0"
        else local cond "`R' < 0"
        qui reg `1' `R' if `cond'
        return scalar delta`j' = _b[`R']
        qui reg `R' `1' if `cond'
        return scalar b`j' = _b[`1']
        qui corr `1' `R' `2' `tau' if `cond', cov
        matrix C = r(C)
        return scalar vX`j' = C[1,1]
        return scalar vR`j' = C[2,2]
        return scalar cXeta`j' = C[1,3]
        return scalar st`j' = C[1,4]/C[2,2]
    }
end

* ---- Panel A: X symmetric, eta homoskedastic (p. 105-107, 110) ----
input X eta
-4 -5
-4  0
-4  5
-3 -5
-3  0
-3  5
 3 -5
 3  0
 3  5
 4 -5
 4  0
 4  5
end
gen R = X + eta
assert R[1] == -9 & R[6] == 2 & R[12] == 9
basupair X eta
scalar A_b = r(b)
scalar A_delta = r(delta)
scalar A_vXf = r(vXf)
scalar A_vRf = r(vRf)
forvalues j = 0/1 {
    scalar A_d`j' = r(delta`j')
    scalar A_b`j' = r(b`j')
    scalar A_ratio`j' = r(vX`j')/r(vR`j')
    scalar A_vX`j' = r(vX`j')
    scalar A_vR`j' = r(vR`j')
    scalar A_st`j' = r(st`j')
    scalar A_cx`j' = r(cXeta`j')
}
di _n "Panel A checks against DMR's printed values (tolerance 0.001):"
di "  full-sample forward slope 1.000 : " %6.3f A_b
di "  full-sample s2X 13.636, s2R 31.818, ratio 0.429 : " %7.3f A_vXf " " %7.3f A_vRf " " %6.3f A_vXf/A_vRf
di "  full-sample reverse slope 0.429  : " %6.3f A_delta
di "  subsample s2X 13.367, s2R 10.7, ratio 1.249 : " %7.3f A_vX0 " " %7.3f A_vR0 " " %6.3f A_ratio0
di "  subsample forward slope 0.651   : " %6.3f A_b0 " " %6.3f A_b1
di "  good and bad reverse slopes 0.813, 0.813 : " %6.3f A_d0 " " %6.3f A_d1
di "  ST ratio with n-1 (-0.436; DMR print -0.363 with n) : " %6.3f A_st0 " " %6.3f A_st1
di "  subsample cov(X, eta) with n (-3.889) : " %6.3f A_cx0*5/6 " " %6.3f A_cx1*5/6
assert abs(A_b - 1.000) < 0.001
assert abs(A_vXf - 13.636) < 0.001 & abs(A_vRf - 31.818) < 0.001
assert abs(A_delta - 0.429) < 0.001
assert abs(A_vX0 - 13.367) < 0.001 & abs(A_vR0 - 10.7) < 0.001 & abs(A_ratio0 - 1.249) < 0.001
assert abs(A_b0 - 0.651) < 0.001 & abs(A_b1 - 0.651) < 0.001
assert abs(A_d0 - 0.813) < 0.001 & abs(A_d1 - 0.813) < 0.001
assert abs(A_d0 - A_d1) < 1e-10
assert abs(A_d0 - (A_b*A_ratio0 + A_st0)) < 1e-10
assert abs(A_cx0*5/6 - (-3.889)) < 0.001

* ---- Panel B: eta variance rising in X (pp. 108, 110) ----
clear
input X eta
-4 -1
-4  0
-4  1
-3 -1
-3  0
-3  1
 3 -5
 3  0
 3  5
 4 -5
 4  0
 4  5
end
gen R = X + eta
assert R[1] == -5 & R[9] == 8 & R[12] == 9
basupair X eta
scalar B_d0 = r(delta0)
scalar B_d1 = r(delta1)
di _n "Panel B checks: good 0.038, bad 1.917 : " %6.3f B_d0 " " %6.3f B_d1
assert abs(B_d0 - 0.038) < 0.001 & abs(B_d1 - 1.917) < 0.001

* ---- Panel C: X left-skewed (pp. 109, 111) ----
clear
input X eta
-6 -5
-6  0
-6  5
-1 -5
-1  0
-1  5
 3 -5
 3  0
 3  5
 4 -5
 4  0
 4  5
end
gen R = X + eta
assert R[1] == -11 & R[6] == 4 & R[12] == 9
basupair X eta
scalar C_d0 = r(delta0)
scalar C_d1 = r(delta1)
di _n "Panel C checks: good 0.281, bad 0.602 : " %6.3f C_d0 " " %6.3f C_d1
assert abs(C_d0 - 0.281) < 0.001 & abs(C_d1 - 0.602) < 0.001

* ---- Mirror of Panel C: X right-skewed (negate X and eta) ----
gen Xm = -X
gen etam = -eta
basupair Xm etam
scalar M_d0 = r(delta0)
scalar M_d1 = r(delta1)
di _n "Mirror identity checks (should be zero to machine precision):"
di "  delta_1(mirror) minus delta_0(C) : " %10.6f (M_d1 - C_d0)
di "  delta_0(mirror) minus delta_1(C) : " %10.6f (M_d0 - C_d1)
assert abs(M_d1 - C_d0) < 1e-10
assert abs(M_d0 - C_d1) < 1e-10

* ---- numbers that travel to the slides ----
foreach s in A_b A_delta A_b0 A_ratio0 A_st0 A_d0 B_d0 B_d1 C_d0 C_d1 M_d0 M_d1 {
    local f`s' : di %6.3f `s'
    local f`s' = strtrim("`f`s''")
}
local fAratiof : di %6.3f A_vXf/A_vRf
local fAratiof = strtrim("`fAratiof'")
* differences computed from the slopes as printed (three decimals), so that
* the table's columns add up as displayed
local fBdiff : di %6.3f round(B_d1, 0.001) - round(B_d0, 0.001)
local fBdiff = strtrim("`fBdiff'")
local fCdiff : di %6.3f round(C_d1, 0.001) - round(C_d0, 0.001)
local fCdiff = strtrim("`fCdiff'")
local fMdiff : di %6.3f round(M_d1, 0.001) - round(M_d0, 0.001)
local fMdiff = strtrim("`fMdiff'")

* Panel A: the components of the full-sample and subsample slopes
file open tab using tables/tab_basu_dmr_a.tex, write replace
file write tab "\begin{tabular}{@{}lccc@{}}" _n
file write tab "\toprule" _n
file write tab " & Full sample & Good, \$R \ge 0\$ & Bad, \$R < 0\$ \\" _n
file write tab "\midrule" _n
file write tab "Full-sample forward slope \$\hat b\$ & `fA_b' & `fA_b' & `fA_b' \\" _n
file write tab "Variance ratio \$s^2_X/s^2_R\$ & `fAratiof' & `fA_ratio0' & `fA_ratio0' \\" _n
file write tab "ST term \$s_{X\tau}/s^2_R\$ & 0.000 & `fA_st0' & `fA_st0' \\" _n
file write tab "Basu slope \$\hat\delta = \hat b \times \text{ratio} + \text{ST}\$ & `fA_delta' & \textcolor{DeepNavy}{\bfseries `fA_d0'} & \textcolor{DeepNavy}{\bfseries `fA_d0'} \\" _n
file write tab "\midrule" _n
file write tab "Subsample forward slope \$\hat b_j\$ (\$\hat\delta_j = \hat b_j \times \text{ratio}_j\$) & `fA_b' & `fA_b0' & `fA_b0' \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_basu_dmr_a.tex
di "Done: tables/tab_basu_dmr_a.tex"

* Panels B, C, and the mirror of C
file open tab using tables/tab_basu_dmr_bc.tex, write replace
file write tab "\begin{tabular}{@{}lccc@{}}" _n
file write tab "\toprule" _n
file write tab "Twelve observations, \$R = X + \eta\$ & Good \$\hat\delta_0\$ & Bad \$\hat\delta_1\$ & \$\hat\beta_1 = \hat\delta_1 - \hat\delta_0\$ \\" _n
file write tab "\midrule" _n
file write tab "A: \$X\$ symmetric, \$\eta\$ homoskedastic & `fA_d0' & `fA_d0' & 0.000 \\" _n
file write tab "B: variance of \$\eta\$ rising in \$X\$ & `fB_d0' & `fB_d1' & \textcolor{Maroon}{\bfseries `fBdiff'} \\" _n
file write tab "C: \$X\$ left-skewed & `fC_d0' & `fC_d1' & \textcolor{Maroon}{\bfseries `fCdiff'} \\" _n
file write tab "Mirror of C: \$X\$ right-skewed, as per-share earnings (DMR fn.~18) & `fM_d0' & `fM_d1' & \textcolor{DeepNavy}{\bfseries `fMdiff'} \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_basu_dmr_bc.tex
di "Done: tables/tab_basu_dmr_bc.tex"


********************************************************************************
* MODULE 4: DMR Table 2, recreated: the placebo return
********************************************************************************
* Week 4, Basu (1997): the placebo return of Dietrich, Muller, and Riedl
* (2007, Table 2, p. 114), recreated
*
* Plots the "bad" and "good" news slopes from three columns of DMR Table 2:
* (1) Basu (1997) Table 2, 1963-1990; (2) DMR's replication on actual returns,
* 1963-1990; (3) the same sample with simulated returns built from the
* full-sample forward regression and another firm's residual, so that the
* return carries no firm-specific asymmetric timeliness. Panels D and E of the
* table decompose each subsample slope as full-sample forward slope + SVR bias
* + ST bias; the script asserts that arithmetic to 0.0015 for columns (2) and
* (3) (column (4) rounds to 0.250 and -0.015 against the printed 0.251 and
* -0.016, DMR extraction flag 4).
* Numbers are typed from the extraction (Table 2, p. 114): no simulation.
* Requirements: base Stata. Output: figures/fig3_basu_table2.png

clear all

input str32 column grp bad good full svrbad stbad svrgood stgood
"Basu (1997) Table 2"        1 0.222 0.061 0.846 . . . .
"DMR, actual returns"        2 0.284 0.025 0.852 -0.055 -0.513 -0.768 -0.059
"DMR, simulated returns"     3 0.493 0.009 0.879 -0.002 -0.384 -0.814 -0.056
end

di _n "DMR Table 2 decomposition checks (full + SVR + ST minus subsample slope, within 0.0015):"
forvalues i = 2/3 {
    di "  " column[`i'] ": bad  " %7.4f (full[`i'] + svrbad[`i'] + stbad[`i'] - bad[`i']) ///
       "   good " %7.4f (full[`i'] + svrgood[`i'] + stgood[`i'] - good[`i'])
    assert abs(full[`i'] + svrbad[`i'] + stbad[`i'] - bad[`i']) < 0.0015
    assert abs(full[`i'] + svrgood[`i'] + stgood[`i'] - good[`i']) < 0.0015
}
di "  Basu Table 2, XE: 0.061 + 0.161 = " %6.3f (0.061 + 0.161) " (bad-news slope 0.222)"
assert abs(0.061 + 0.161 - bad[1]) < 0.0015

* Like with like: the asymmetric timeliness coefficient beta1 is the bad-news
* slope minus the good-news slope in each column, which is the quantity to set
* beside Basu's own beta1 (Basu Table 1 Panel A: 0.216; Basu Table 2, XE: 0.161).
* Slope ratios (Basu's 4.66 and 3.64) are a different quantity and are not
* compared with a simulated coefficient.
di _n "Asymmetric timeliness coefficient beta1 = bad-news slope minus good-news slope:"
forvalues i = 1/3 {
    di "  " column[`i'] ": " %5.3f bad[`i'] " - " %5.3f good[`i'] " = " %6.3f (bad[`i'] - good[`i'])
}
di "  Basu Table 1 Panel A (printed beta1): 0.216;  Basu Table 2, XE (printed beta1): 0.161"
assert abs((bad[1] - good[1]) - 0.161) < 0.0015
assert (bad[3] - good[3]) > 0.216

gen xb = grp - 0.19
gen xg = grp + 0.19
gen lb = string(bad, "%5.3f")
gen lg = string(good, "%5.3f")
gen yb = bad + 0.02
gen yg = good + 0.02

twoway (bar bad xb, barwidth(0.36) color(maroon))                                        ///
       (bar good xg, barwidth(0.36) color(navy))                                          ///
       (scatter yb xb, msymbol(none) mlabel(lb) mlabposition(12) mlabcolor(maroon) mlabsize(medsmall)) ///
       (scatter yg xg, msymbol(none) mlabel(lg) mlabposition(12) mlabcolor(navy) mlabsize(medsmall)), ///
    xlabel(1 "Basu (1997) Table 2" 2 "DMR replication, actual returns" 3 "DMR, returns rebuilt with no asymmetry", labsize(small)) ///
    xscale(range(0.5 3.5)) yscale(range(0 0.58))                                          ///
    ylabel(0(0.1)0.5, angle(horizontal) format(%3.1f))                                    ///
    ytitle("Slope of scaled earnings on returns") xtitle("")                              ///
    title("Returns built with no asymmetry reproduce the gap")                            ///
    subtitle("Recreated from Dietrich, Muller, and Riedl (2007), Table 2, columns (1) to (3), 1963-1990") ///
    legend(order(1 "bad news, R < 0" 2 "good news, R {&ge} 0") ring(0) position(11)       ///
           cols(1) size(medsmall) region(lstyle(none)))                                   ///
    scheme(s1color) xsize(9) ysize(5)
graph export figures/fig3_basu_table2.png, replace width(2000)
di "Done: figures/fig3_basu_table2.png"


********************************************************************************
* MODULE 5: The Basu coefficient under the null across eight designs
********************************************************************************
* Week 4, Basu (1997): the sampling distribution of the asymmetric timeliness
* coefficient when the true asymmetry is zero, across eight designs
*
* Demonstrates that under the null of no conservatism the Basu interaction
* coefficient beta1 = delta_1 - delta_0 is centered at zero only when
* Dietrich, Muller, and Riedl's (2007, pp. 104-105) sufficient conditions hold,
* and otherwise takes a sign and a magnitude set by the joint distribution of
* earnings and non-earnings information and by the location of the split.
*
* DGP (DMR eq. 1.1, beta = lambda = 1, cov(X, eta) = 0): R = X + eta, moments
* calibrated to DMR Table 1: sd(X) = 0.17, sd(eta) = 0.497 (var(R) = 0.276).
* Eight designs, 1,000 replications of N = 2,000 each, split at R = 0:
*   A  X, eta symmetric normal, mean zero (split at the mean; conditions hold)
*   B  eta right-skewed: standardized chi-square(2), mean 0 (returns' right tail)
*   C  X left-skewed: minus a standardized chi-square(2) (earnings' left tail)
*   D  X right-skewed (per-share earnings; DMR Table 3 Panel D, fn. 18)
*   E  eta heteroskedastic: sd(eta | X) = 0.497 sqrt(2 Phi(X/0.17)), rising in X
*   F  X, eta normal with mean X = 0.063, mean eta = 0.101, so the mean of R is
*      0.164 as in DMR Table 1 and the split at zero is off the mean; under joint
*      normality E[X | R] is linear, so no bias is expected whatever the split
*   G  X, eta symmetric Laplace with the same means, split at zero: symmetric
*      but non-normal data split away from the mean
*   H  DMR Table 1's shape: X left-skewed with mean 0.063, eta right-skewed with
*      mean 0.101, split at zero
* Each replication records the two subsample slopes, beta1 and its t-statistic
* from the stacked regression, and the two subsample R2s.
* Checks: designs A and F mean beta1 within 0.01 of zero (about three Monte
* Carlo standard errors); design D mean beta1 equals minus design C mean within
* 0.01 (mirror).
* Requirements: base Stata. Runtime about two minutes.
* Outputs: figures/fig2_basu_designs.png, tables/tab_basu_montecarlo.tex

clear all
set seed 1997

local reps = 1000
local N    = 2000

tempname sims
postfile `sims' str1 design rep d0 d1 phi1 t1 r2b r2g using mc_basu_results, replace
foreach d in A B C D E F G H {
    forvalues r = 1/`reps' {
        quietly {
            clear
            set obs `N'
            if "`d'" == "A" {
                gen X   = 0.17*rnormal()
                gen eta = 0.497*rnormal()
            }
            if "`d'" == "B" {
                gen X   = 0.17*rnormal()
                gen eta = 0.497*(rchi2(2) - 2)/2
            }
            if "`d'" == "C" {
                gen X   = -0.17*(rchi2(2) - 2)/2
                gen eta = 0.497*rnormal()
            }
            if "`d'" == "D" {
                gen X   = 0.17*(rchi2(2) - 2)/2
                gen eta = 0.497*rnormal()
            }
            if "`d'" == "E" {
                gen X   = 0.17*rnormal()
                gen eta = 0.497*rnormal()*sqrt(2*normal(X/0.17))
            }
            if "`d'" == "F" {
                gen X   = 0.063 + 0.17*rnormal()
                gen eta = 0.101 + 0.497*rnormal()
            }
            if "`d'" == "G" {
                gen X   = 0.063 + 0.17*(-ln(runiform()) + ln(runiform()))/sqrt(2)
                gen eta = 0.101 + 0.497*(-ln(runiform()) + ln(runiform()))/sqrt(2)
            }
            if "`d'" == "H" {
                gen X   = 0.063 - 0.17*(rchi2(2) - 2)/2
                gen eta = 0.101 + 0.497*(rchi2(2) - 2)/2
            }
            gen R   = X + eta
            gen DR  = (R < 0)
            gen RDR = R*DR
            reg X DR R RDR
            local phi1 = _b[RDR]
            local t1   = _b[RDR]/_se[RDR]
            reg X R if R >= 0
            local d0  = _b[R]
            local r2g = e(r2)
            reg X R if R < 0
            local d1  = _b[R]
            local r2b = e(r2)
            post `sims' ("`d'") (`r') (`d0') (`d1') (`phi1') (`t1') (`r2b') (`r2g')
        }
    }
}
postclose `sims'

use mc_basu_results, clear
gen rej    = abs(t1) > 1.96
gen r2diff = r2b - r2g
gen lo = .
gen hi = .
foreach d in A B C D E F G H {
    _pctile phi1 if design == "`d'", p(2.5 97.5)
    replace lo = r(r1) if design == "`d'"
    replace hi = r(r2) if design == "`d'"
}
collapse (mean) d0 d1 phi1 t1 rej r2diff lo hi (sd) sd_phi1 = phi1, by(design)
di _n "Monte Carlo summary by design (`reps' replications of N = `N', true beta1 = 0):"
format d0 d1 phi1 t1 rej r2diff sd_phi1 lo hi %8.3f
list design d0 d1 phi1 sd_phi1 lo hi t1 rej r2diff, clean noobs

* ---- checks ----
qui su phi1 if design == "A"
scalar mA = r(mean)
qui su sd_phi1 if design == "A"
scalar seA = r(mean)/sqrt(`reps')
qui su phi1 if design == "F"
scalar mF = r(mean)
qui su phi1 if design == "C"
scalar mC = r(mean)
qui su phi1 if design == "D"
scalar mD = r(mean)
di _n "Monte Carlo checks:"
di "  design A mean beta1 (should be within 0.01 of zero; MC s.e. " %6.4f seA "): " %8.4f mA
di "  design F mean beta1 (normal, split off the mean; within 0.01 of zero)  : " %8.4f mF
di "  design D mean beta1 plus design C mean beta1 (mirror; within 0.01)     : " %8.4f (mD + mC)
assert abs(mA) < 0.01
assert abs(mF) < 0.01
assert abs(mD + mC) < 0.01

* ---- table fragment ----
gen ypos = .
replace ypos = 8 if design == "A"
replace ypos = 7 if design == "B"
replace ypos = 6 if design == "C"
replace ypos = 5 if design == "D"
replace ypos = 4 if design == "E"
replace ypos = 3 if design == "F"
replace ypos = 2 if design == "G"
replace ypos = 1 if design == "H"
sort ypos
gen rejlab = string(100*rej, "%3.0f") + "% reject"

file open tab using tables/tab_basu_montecarlo.tex, write replace
file write tab "\begin{tabular}{@{}lcccc@{}}" _n
file write tab "\toprule" _n
file write tab "Design (true \$\beta_1 = 0\$) & Mean \$\hat\beta_1\$ & 2.5th, 97.5th pct. & Reject at 5\% & Mean \$R^2_{\mathrm{bad}} - R^2_{\mathrm{good}}\$ \\" _n
file write tab "\midrule" _n
gsort -ypos
local texlab1 "A: normal, split at the mean"
local texlab2 "B: \(\eta\) right-skewed"
local texlab3 "C: \(X\) left-skewed"
local texlab4 "D: \(X\) right-skewed"
local texlab5 "E: \(\mathrm{var}(\eta)\) rising in \(X\)"
local texlab6 "F: normal, split at 0, mean \(R\) = 0.16"
local texlab7 "G: heavy-tailed, split at 0, mean \(R\) = 0.16"
local texlab8 "H: DMR Table 1 shape, split at 0"
forvalues i = 1/8 {
    local m  : di %6.3f phi1[`i']
    local m  = strtrim("`m'")
    if abs(phi1[`i']) < 0.0005 local m "0.000"
    local l  : di %6.3f lo[`i']
    local l  = strtrim("`l'")
    local h  : di %6.3f hi[`i']
    local h  = strtrim("`h'")
    local rj : di %3.0f 100*rej[`i']
    local rj = strtrim("`rj'")
    local r2 : di %6.3f r2diff[`i']
    local r2 = strtrim("`r2'")
    if `i' == 1 | `i' == 6 local cell "`m'"
    else if phi1[`i'] > 0 local cell "\textcolor{Maroon}{\bfseries `m'}"
    else local cell "\textcolor{DeepNavy}{\bfseries `m'}"
    file write tab "`texlab`i'' & `cell' & `l', `h' & `rj'\% & `r2' \\" _n
}
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_basu_montecarlo.tex
di "Done: tables/tab_basu_montecarlo.tex"

* ---- figure: mean and 2.5th to 97.5th percentile band by design ----
gen xlab = cond(hi > 0, hi + 0.03, 0.03)
twoway (rcap lo hi ypos, horizontal lcolor(gs8) lwidth(medthick) msize(medium))         ///
       (scatter ypos phi1 if phi1 >= 0.01, mcolor(maroon) msymbol(circle) msize(large)) ///
       (scatter ypos phi1 if phi1 <= -0.01, mcolor(navy) msymbol(circle) msize(large)) ///
       (scatter ypos phi1 if abs(phi1) < 0.01, mcolor(gs6) msymbol(circle) msize(large)) ///
       (scatter ypos xlab, msymbol(none) mlabel(rejlab) mlabcolor(gs6) mlabsize(small) mlabposition(3)), ///
    xline(0, lcolor(gs8) lpattern(dash))                                                 ///
    ylabel(8 "A: normal, split at the mean" 7 "B: {&eta} right-skewed"                   ///
           6 "C: X left-skewed" 5 "D: X right-skewed" 4 "E: var({&eta}) rising in X"       ///
           3 "F: normal, split at 0, mean R = 0.16" 2 "G: heavy-tailed, same split"        ///
           1 "H: DMR Table 1 shape, split at 0",                                              ///
           angle(horizontal) labsize(small) nogrid)                                      ///
    ytitle("") xtitle("Estimated asymmetric timeliness coefficient, true value 0")     ///
    xlabel(-0.4(0.2)1.0, format(%3.1f)) xscale(range(-0.55 1.2))                         ///
    yscale(range(0.35 8.5))                                                              ///
    text(0.55 0.03 "positive: read as conservative", color(maroon) size(small) placement(e)) ///
    text(0.55 -0.03 "negative: read as aggressive", color(navy) size(small) placement(w)) ///
    title("Under the null the coefficient takes either sign")                           ///
    subtitle("1,000 samples of N = 2,000 (seed 1997); bars: 2.5th to 97.5th percentiles") ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig2_basu_designs.png, replace width(2000)
di "Done: figures/fig2_basu_designs.png"
erase mc_basu_results.dta

log close
