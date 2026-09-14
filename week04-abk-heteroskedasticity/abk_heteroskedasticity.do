********************************************************************************
* ACCTG 597A -- Week 4: Heteroskedasticity of Known and Unknown Form
* Aboody, Barth, and Kasznik (1999), "Revaluations of Fixed Assets and Future
* Firm Performance: Evidence from the UK," Journal of Accounting and Economics
* 26(1-3): 149-178. Equation (1), p. 158.
*
* Self-contained simulation, base Stata only (version 17 or later). Run from
* this folder:  do abk_heteroskedasticity.do
* Figures are written to figures/ and booktabs table fragments to tables/.
* Runtime about three minutes (Module 3 runs 10,000 samples of N = 740).
*
* Modules. In each, the true slope beta = 0.07 is known by construction, so
* the property under study is visible directly. Seed 1999 throughout, so the
* numbers match the lecture deck.
*   1. Identities in one sample. Levels errors with s.d. 0.065 MVE. The
*      conventional OLS variance is wrong by 1 + Cov(x~^2, sigma^2)/(mean x~^2
*      mean sigma^2) = 6.16; White's estimator equals the sandwich with squared
*      residuals and Stata's vce(robust) to eight decimals; deflating every
*      variable by MVE equals regress [aweight = 1/MVE^2] exactly.
*   2. Residual fans. The same draws in levels (fan), deflated at the right
*      exponent (flat band), and deflated at the wrong exponent (reverse fan).
*   3. Monte Carlo. Ten designs (error s.d. proportional to MVE^(gamma/2) for
*      gamma = 0 to 4, plus a regressor-dependent variance no deflator removes);
*      size of conventional, HC1, and HC3 t-tests in levels and after deflation;
*      relative efficiency of the deflated slope (6.5 at gamma = 0, 0.21 at
*      gamma = 2).
********************************************************************************

version 17
clear all
set more off
capture log close
log using abk_heteroskedasticity.log, replace text
capture mkdir figures
capture mkdir tables


********************************************************************************
* MODULE 1: Identities in one sample: the covariance ratio, White = sandwich with squared residuals, deflation = WLS
********************************************************************************
********************************************************************************
* ACCTG 597A -- Week 4: Heteroskedasticity of Known and Unknown Form
* Aboody, Barth, and Kasznik (1999), JAE 26: 149-178, equation (1)
*
* Script 1: identities in one sample.
* Demonstrates (i) that the conventional OLS variance is wrong by the factor
* 1 + Cov(x~^2, sigma^2)/(mean x~^2 * mean sigma^2), computed here with the
* error variances known by construction (x~ is the regressor after partialling
* out the other regressors); (ii) that White's (1980) estimator is the sandwich
* formula with squared residuals in place of the unknown error variances and
* equals Stata's vce(robust) up to the finite-sample factor n/(n-k); (iii) that
* deflating every variable, including the constant, by scale z is weighted
* least squares with analytic weight 1/z^2.
*
* DGP, calibrated to ABK Table 2 (p. 161). One firm-year per observation.
*   z = beginning market value of equity (MVE), lognormal, log z ~ N(4, 1).
*   r = REV/MVE, the revaluation ratio: r = rgamma(0.35, 0.314) (z/z_med)^(-0.25),
*       mean about 0.11, median about 0.03, s.d. about 0.19, declining in size.
*   x = REV = z r (pounds).
*   y = change in operating income, y = 0.02 z + 0.07 x + u, u = 0.065 z e,
*       e ~ N(0, 1): sigma_i = 0.065 z_i, heteroskedasticity of known form with
*       exponent gamma = 2, so the deflated regression y/z on x/z is
*       homoskedastic with error s.d. 0.065 (ABK: s.d. of dOPINC/MVE = 0.07).
*   True slope beta = 0.07 (ABK Table 3, one-year horizon). N = 740.
* The levels regression is y on x, z, and a constant (the constant is zero in
* truth; the levels coefficient on z is the deflated intercept). The deflated
* regression is y/z on x/z, 1/z, and a constant.
*
* Requirements: base Stata (version 17 or later).
* Output: tables/tab_abk_heteroskedasticity_identity.tex
********************************************************************************
clear all
set seed 1999

local n = 740
set obs `n'
gen z = exp(4 + rnormal())
quietly summarize z, detail
scalar zmed = r(p50)
gen r = rgamma(0.35, 0.314)*(z/zmed)^(-0.25)
gen x = z*r
gen s2 = (0.065*z)^2                    // sigma_i^2, known by construction
gen u = sqrt(s2)*rnormal()
gen y = 0.02*z + 0.07*x + u
quietly summarize r, detail
di _n "Calibration: r = REV/MVE mean " %5.3f r(mean) ", median " %5.3f r(p50) ", s.d. " %5.3f r(sd)
quietly summarize z, detail
di "             z = MVE median " %6.1f r(p50) ", p90/p10 " %5.1f r(p90)/r(p10)

* ---- (i) OLS in levels: conventional and White standard errors ----
quietly regress y x z
scalar b_ols   = _b[x]
scalar se_conv = _se[x]
scalar k       = e(df_m) + 1
predict uhat, resid
quietly regress y x z, vce(robust)
scalar se_white = _se[x]
quietly regress y x z, vce(hc3)
scalar se_hc3 = _se[x]

* Hand-computed White (HC0) for the slope on x: sum(xt^2 uhat^2)/(sum xt^2)^2,
* xt the part of x orthogonal to z and the constant (Frisch-Waugh-Lovell).
quietly regress x z
predict xt, resid
quietly generate double num = xt^2*uhat^2
quietly summarize num
scalar S_num = r(sum)
quietly generate double den = xt^2
quietly summarize den
scalar Sxx = r(sum)
scalar v_hc0  = S_num/Sxx^2
scalar se_hc1 = sqrt(v_hc0*`n'/(`n' - k))

* True sampling variance of the slope, using the known sigma_i^2:
quietly generate double tnum = xt^2*s2
quietly summarize tnum
scalar v_true  = r(sum)/Sxx^2
scalar se_true = sqrt(v_true)
* Conventional target sigma_bar^2/Sxx and the covariance identity:
quietly summarize s2
scalar sbar2 = r(mean)
scalar v_convtarget = sbar2/Sxx
quietly correlate den s2, covariance
scalar cov_x2s2 = r(cov_12)*(`n' - 1)/`n'      // divisor n
scalar mx2 = Sxx/`n'
scalar ratio_identity = 1 + cov_x2s2/(mx2*sbar2)

* ---- (ii) deflation by z equals WLS with weight 1/z^2 ----
gen yd   = y/z
gen xd   = x/z
gen invz = 1/z
quietly regress yd xd invz
scalar b_defl   = _b[xd]
scalar se_defl  = _se[xd]
scalar a_defl   = _b[_cons]
quietly regress yd xd invz, vce(robust)
scalar se_deflw = _se[xd]
gen w = 1/z^2
quietly regress y x z [aweight = w]
scalar b_wls  = _b[x]
scalar se_wls = _se[x]
scalar a_wls  = _b[z]

* ---- identity checks: print them for the log, then assert them ----
di _n "Identity checks (should be zero to machine precision):"
di "  Stata vce(robust) minus hand HC1 (HC0 * n/(n-k)) : " %12.8f (se_white - se_hc1)
di "  deflated OLS slope minus WLS [aw=1/z^2] slope   : " %12.8f (b_defl - b_wls)
di "  deflated OLS s.e. minus WLS [aw=1/z^2] s.e.     : " %12.8f (se_defl - se_wls)
di "  deflated intercept minus WLS coefficient on z   : " %12.8f (a_defl - a_wls)
di "  v_true/v_convtarget minus covariance identity   : " %12.8f (v_true/v_convtarget - ratio_identity)
assert abs(se_white - se_hc1) < 1e-6
assert abs(b_defl - b_wls) < 1e-6
assert abs(se_defl - se_wls) < 1e-6
assert abs(a_defl - a_wls) < 1e-6
assert abs(v_true/v_convtarget - ratio_identity) < 1e-6

di _n "Numbers for the slide:"
di "  levels OLS slope           : " %8.4f b_ols
di "  conventional s.e.          : " %8.4f se_conv
di "  White (HC1) s.e.           : " %8.4f se_white
di "  White (HC3) s.e.           : " %8.4f se_hc3
di "  true s.e. (known sigma_i)  : " %8.4f se_true
di "  ratio true/conventional var: " %8.3f ratio_identity
di "  deflated slope             : " %8.4f b_defl
di "  deflated conventional s.e. : " %8.4f se_defl
di "  deflated White s.e.        : " %8.4f se_deflw

* ---- numbers that travel to the slide ----
foreach s in b_ols se_conv se_white se_hc3 se_true b_defl se_defl se_deflw {
    local f`s' : di %5.3f `s'
    local f`s' = strtrim("`f`s''")
}
local t_conv  : di %3.1f b_ols/se_conv
local t_conv  = strtrim("`t_conv'")
local t_white : di %3.1f b_ols/se_white
local t_white = strtrim("`t_white'")
local t_hc3   : di %3.1f b_ols/se_hc3
local t_hc3   = strtrim("`t_hc3'")
local t_defl  : di %3.1f b_defl/se_defl
local t_defl  = strtrim("`t_defl'")
local t_deflw : di %3.1f b_defl/se_deflw
local t_deflw = strtrim("`t_deflw'")
local fratio : di %3.1f ratio_identity
local fratio = strtrim("`fratio'")

* ---- booktabs fragment (escape math dollars as \$) ----
capture mkdir tables
file open tab using tables/tab_abk_heteroskedasticity_identity.tex, write replace
file write tab "\begin{tabular}{@{}llcc@{}}" _n
file write tab "\toprule" _n
file write tab "Regression & Standard error & \$\hat\beta\$ & s.e. (\$t\$) \\" _n
file write tab "\midrule" _n
file write tab "Levels & conventional \$s^2/S_{xx}\$ & `fb_ols' & \textcolor{Maroon}{\bfseries `fse_conv' (`t_conv')} \\" _n
file write tab "Levels & White (1980), HC1 & `fb_ols' & `fse_white' (`t_white') \\" _n
file write tab "Levels & White, HC3 & `fb_ols' & `fse_hc3' (`t_hc3') \\" _n
file write tab "Levels & true, \$\sigma_i\$ known & `fb_ols' & \textcolor{DeepNavy}{\bfseries `fse_true'} \\" _n
file write tab "\midrule" _n
file write tab "Deflated by \$z\$ & conventional & `fb_defl' & \textcolor{DeepNavy}{\bfseries `fse_defl' (`t_defl')} \\" _n
file write tab "Deflated by \$z\$ & White (1980), HC1 & `fb_defl' & `fse_deflw' (`t_deflw') \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_abk_heteroskedasticity_identity.tex
di "Done: tables/tab_abk_heteroskedasticity_identity.tex"

********************************************************************************
* MODULE 2: The picture: residual fans in levels and after deflation, at the right and the wrong exponent
********************************************************************************
********************************************************************************
* ACCTG 597A -- Week 4: Heteroskedasticity of Known and Unknown Form
* Aboody, Barth, and Kasznik (1999), JAE 26: 149-178, equation (1)
*
* Script 2: the picture. Residuals against scale for the levels regression
* (a fan that widens with z when sigma_i = 0.065 z) and for the same data after
* every variable is deflated by z (a band of constant width). A third panel
* deflates data whose error standard deviation is proportional to the square
* root of z (gamma = 1): deflation by z now over-corrects and the band narrows
* with z, so the deflated regression is still heteroskedastic.
*
* DGP as in script 1 (seed 1999, N = 740; calibrated to ABK Table 2):
* log z ~ N(4, 1); r = rgamma(0.35, 0.314) (z/z_med)^(-0.25); x = z r;
* y = 0.02 z + 0.07 x + u; u = sigma_gamma z^(gamma/2) e, e ~ N(0, 1), with
* sigma_gamma = 0.065 z_med^((2 - gamma)/2), so the error s.d. equals 0.065 z
* at the median firm for every gamma. Levels regression: y on x, z, constant.
* Deflated regression: y/z on x/z, 1/z, constant.
*
* Requirements: base Stata (version 17 or later).
* Output: figures/fig1_abk_heteroskedasticity_fan.png
********************************************************************************
clear all
set seed 1999

local n = 740
set obs `n'
gen z = exp(4 + rnormal())
quietly summarize z, detail
scalar zmed = r(p50)
gen r = rgamma(0.35, 0.314)*(z/zmed)^(-0.25)
gen x = z*r
gen e = rnormal()
gen lz = log(z)

* gamma = 2: sigma_i = 0.065 z
gen u2 = 0.065*z*e
gen y2 = 0.02*z + 0.07*x + u2
quietly regress y2 x z
predict res_lev, resid
gen y2d = y2/z
gen xd  = x/z
gen invz = 1/z
quietly regress y2d xd invz
predict res_d2, resid

* gamma = 1: sigma_i = sigma_1 sqrt(z), then deflated by z
scalar sig1 = 0.065*zmed^(0.5)
gen u1 = sig1*sqrt(z)*e
gen y1 = 0.02*z + 0.07*x + u1
gen y1d = y1/z
quietly regress y1d xd invz
predict res_d1, resid

di _n "Checks: the error s.d. at the median firm is 0.065 z_med for both gammas"
di "  gamma = 2: 0.065*zmed       = " %10.4f 0.065*zmed
di "  gamma = 1: sig1*sqrt(zmed)  = " %10.4f sig1*sqrt(zmed)
assert abs(0.065*zmed - sig1*sqrt(zmed)) < 1e-6

* Sample s.d. of deflated residuals in the bottom and top scale terciles
xtile ter = z, nq(3)
foreach v in res_d2 res_d1 {
    quietly summarize `v' if ter == 1
    scalar sd1_`v' = r(sd)
    quietly summarize `v' if ter == 3
    scalar sd3_`v' = r(sd)
    di "  `v': s.d. bottom tercile " %7.4f sd1_`v' "  top tercile " %7.4f sd3_`v' ///
       "  ratio top/bottom " %5.2f sd3_`v'/sd1_`v'
}
local r2 : di %3.1f sd3_res_d2/sd1_res_d2
local r2 = strtrim("`r2'")
local r1 : di %3.1f sd3_res_d1/sd1_res_d1
local r1 = strtrim("`r1'")

twoway (scatter res_lev lz, mcolor(navy%45) msize(vsmall) msymbol(circle)), ///
    yline(0, lcolor(gs8) lpattern(dash)) ///
    title("Levels: residual spread grows with scale", size(medium) color(black)) ///
    subtitle("{&sigma}{sub:i} = 0.065 MVE{sub:i} (known form)", size(small)) ///
    xtitle("log market value of equity") ytitle("Residual (levels)") ///
    ylabel(, angle(horizontal) format(%4.0f)) xlabel(1(1)7) ///
    legend(off) scheme(s1color) name(g1, replace) nodraw
twoway (scatter res_d2 lz, mcolor(navy%45) msize(vsmall) msymbol(circle)), ///
    yline(0, lcolor(gs8) lpattern(dash)) ///
    title("Deflated by MVE: constant spread", size(medium) color(black)) ///
    subtitle("s.d. ratio, top to bottom size tercile: `r2'", size(small)) ///
    xtitle("log market value of equity") ytitle("Deflated residual") ///
    ylabel(-0.2(0.1)0.2, angle(horizontal) format(%3.1f)) xlabel(1(1)7) ///
    legend(off) scheme(s1color) name(g2, replace) nodraw
twoway (scatter res_d1 lz, mcolor(maroon%45) msize(vsmall) msymbol(circle)), ///
    yline(0, lcolor(gs8) lpattern(dash)) ///
    title("Deflated by MVE when {&sigma}{sub:i} {&prop} MVE{sup:1/2}", size(medium) color(black)) ///
    subtitle("over-corrected; s.d. ratio, top to bottom size tercile: `r1'", size(small)) ///
    xtitle("log market value of equity") ytitle("Deflated residual") ///
    ylabel(-0.4(0.2)0.4, angle(horizontal) format(%3.1f)) xlabel(1(1)7) ///
    legend(off) scheme(s1color) name(g3, replace) nodraw
graph combine g1 g2 g3, rows(1) xsize(13) ysize(4.6) scheme(s1color) ///
    title("Deflation removes heteroskedasticity only when its form matches the deflator", size(medsmall) color(black)) ///
    subtitle("Simulated data (seed 1999), N = 740, calibrated to ABK Table 2; the same draws in all three panels", size(small) color(gs6))
capture mkdir figures
graph export figures/fig1_abk_heteroskedasticity_fan.png, replace width(2000)
di "Done: figures/fig1_abk_heteroskedasticity_fan.png"

********************************************************************************
* MODULE 3: Monte Carlo over the form of heteroskedasticity: size, s.e. accuracy, and efficiency of levels versus deflated slopes
********************************************************************************
********************************************************************************
* ACCTG 597A -- Week 4: Heteroskedasticity of Known and Unknown Form
* Aboody, Barth, and Kasznik (1999), JAE 26: 149-178, equation (1)
*
* Script 3: Monte Carlo over the form of heteroskedasticity.
* For each design, 1,000 samples of N = 740 are drawn, the levels regression
* (y on x, z, constant) and the deflated regression (y/z on x/z, 1/z, constant)
* are estimated, and t-tests of the true slope beta = 0.07 are recorded with
* conventional, White (1980) HC1, and HC3 standard errors. Reported: the
* empirical size of a nominal 5 percent test, the ratio of the mean reported
* s.e. to the sampling s.d. of the slope, and the efficiency of the deflated
* slope relative to the levels slope (ratio of sampling s.d.).
*
* Designs. Known-form family: sigma_i = sigma_gamma z_i^(gamma/2) for gamma in
* {0, 0.5, ..., 4}; gamma = 2 is the form that deflation by z assumes and
* gamma = 0 is homoskedastic levels. sigma_gamma = 0.065 z_med^((2 - gamma)/2)
* equalizes the error s.d. at the median firm across gamma. Unknown form:
* sigma_i = 0.065 z_i (1 + 2 D_i), D_i = 1 if r_i exceeds its sample median,
* so the error s.d. is three times larger for above-median revaluers: the
* variance depends on the deflated regressor itself and no scale deflator
* removes it.
* DGP otherwise as in script 1 (calibrated to ABK Table 2): log z ~ N(4, 1),
* r = rgamma(0.35, 0.314) (z/z_med)^(-0.25), x = z r, y = 0.02 z + 0.07 x + u.
* Seed 1999.
* Checks. A correctly specified conventional test (levels at gamma = 0;
* deflated at gamma = 2) has empirical size within 0.02 of 0.05 (binomial s.e.
* 0.007 with 1,000 replications); the deflated White test has size below 0.08
* in every design; both slopes are unbiased for 0.07 in every design (Monte
* Carlo mean within four Monte Carlo standard errors of 0.07).
*
* Requirements: base Stata (version 17 or later). Runtime a few minutes.
* Output: tables/tab_abk_heteroskedasticity_montecarlo.tex,
*         figures/fig2_abk_heteroskedasticity_size.png,
*         figures/fig3_abk_heteroskedasticity_efficiency.png
********************************************************************************
clear all
set seed 1999

local n    = 740
local reps = 1000
local beta = 0.07

tempname sims
tempfile results
postfile `sims' design gamma b_lev se_lev_c se_lev_w se_lev_h3 ///
                b_def se_def_c se_def_w se_def_h3 using `results', replace

* Designs 1 to 9: gamma grid; design 10: unknown form (regressor-dependent).
local gammas 0 0.5 1 1.5 2 2.5 3 3.5 4
local d = 0
foreach g of local gammas {
    local ++d
    forvalues rep = 1/`reps' {
        quietly {
            clear
            set obs `n'
            gen z = exp(4 + rnormal())
            summarize z, detail
            scalar zmed = r(p50)
            gen r = rgamma(0.35, 0.314)*(z/zmed)^(-0.25)
            gen x = z*r
            scalar sg = 0.065*zmed^((2 - `g')/2)
            gen u = sg*z^(`g'/2)*rnormal()
            gen y = 0.02*z + `beta'*x + u
            regress y x z
            scalar b1 = _b[x]
            scalar s1 = _se[x]
            regress y x z, vce(robust)
            scalar s2 = _se[x]
            regress y x z, vce(hc3)
            scalar s3 = _se[x]
            gen yd = y/z
            gen xd = x/z
            gen invz = 1/z
            regress yd xd invz
            scalar b2 = _b[xd]
            scalar s4 = _se[xd]
            regress yd xd invz, vce(robust)
            scalar s5 = _se[xd]
            regress yd xd invz, vce(hc3)
            scalar s6 = _se[xd]
            post `sims' (`d') (`g') (b1) (s1) (s2) (s3) (b2) (s4) (s5) (s6)
        }
    }
    di "design `d' (gamma = `g') done"
}
local ++d
forvalues rep = 1/`reps' {
    quietly {
        clear
        set obs `n'
        gen z = exp(4 + rnormal())
        summarize z, detail
        scalar zmed = r(p50)
        gen r = rgamma(0.35, 0.314)*(z/zmed)^(-0.25)
        gen x = z*r
        summarize r, detail
        gen D = r > r(p50)
        gen u = 0.065*z*(1 + 2*D)*rnormal()
        gen y = 0.02*z + `beta'*x + u
        regress y x z
        scalar b1 = _b[x]
        scalar s1 = _se[x]
        regress y x z, vce(robust)
        scalar s2 = _se[x]
        regress y x z, vce(hc3)
        scalar s3 = _se[x]
        gen yd = y/z
        gen xd = x/z
        gen invz = 1/z
        regress yd xd invz
        scalar b2 = _b[xd]
        scalar s4 = _se[xd]
        regress yd xd invz, vce(robust)
        scalar s5 = _se[xd]
        regress yd xd invz, vce(hc3)
        scalar s6 = _se[xd]
        post `sims' (`d') (.) (b1) (s1) (s2) (s3) (b2) (s4) (s5) (s6)
    }
}
di "design `d' (unknown form) done"
postclose `sims'

use `results', clear
scalar tcrit = invttail(`n' - 3, 0.025)
foreach e in lev_c lev_w lev_h3 def_c def_w def_h3 {
    local b = cond(substr("`e'", 1, 3) == "lev", "b_lev", "b_def")
    gen rej_`e' = abs((`b' - `beta')/se_`e') > tcrit
}
collapse (mean) rej_* m_se_lev_c = se_lev_c m_se_lev_w = se_lev_w m_se_lev_h3 = se_lev_h3 ///
                m_se_def_c = se_def_c m_se_def_w = se_def_w m_se_def_h3 = se_def_h3 ///
                mb_lev = b_lev mb_def = b_def ///
         (sd)   sd_lev = b_lev sd_def = b_def, by(design gamma)
foreach e in lev_c lev_w lev_h3 {
    gen ratio_`e' = m_se_`e'/sd_lev
}
foreach e in def_c def_w def_h3 {
    gen ratio_`e' = m_se_`e'/sd_def
}
gen releff = sd_def/sd_lev
format rej_* ratio_* releff %6.3f
format mb_* %7.4f
list design gamma rej_lev_c rej_lev_w rej_lev_h3 rej_def_c rej_def_w rej_def_h3, noobs sep(0)
list design gamma ratio_lev_c ratio_lev_w ratio_lev_h3 ratio_def_c ratio_def_w ratio_def_h3, noobs sep(0)
list design gamma releff mb_lev mb_def sd_lev sd_def, noobs sep(0)

* ---- checks against the analytic expectations ----
di _n "Checks (a correctly specified conventional test has size within 0.02 of 0.05):"
quietly summarize rej_lev_c if gamma == 0
di "  levels, conventional, gamma = 0   : " %6.3f r(mean)
assert abs(r(mean) - 0.05) < 0.02
quietly summarize rej_def_c if gamma == 2
di "  deflated, conventional, gamma = 2 : " %6.3f r(mean)
assert abs(r(mean) - 0.05) < 0.02
quietly summarize rej_def_w
di "  deflated, White, max over designs  : " %6.3f r(max)
assert r(max) < 0.08
quietly summarize rej_lev_w
di "  levels, White HC1, max over designs (finite-sample leverage; printed, not asserted): " %6.3f r(max)
quietly summarize rej_lev_h3
di "  levels, HC3, max over designs (printed, not asserted): " %6.3f r(max)
di "Checks (both slopes unbiased for 0.07 in every design: the Monte Carlo mean"
di "        lies within four Monte Carlo standard errors, 4 sd/sqrt(1000), of 0.07):"
gen dev_lev = abs(mb_lev - `beta')/(sd_lev/sqrt(`reps'))
gen dev_def = abs(mb_def - `beta')/(sd_def/sqrt(`reps'))
format dev_* %6.2f
list design gamma mb_lev dev_lev mb_def dev_def, noobs sep(0)
quietly summarize dev_lev
di "  levels slope, max deviation in MC s.e. units   : " %6.2f r(max)
assert r(max) < 4
quietly summarize dev_def
di "  deflated slope, max deviation in MC s.e. units : " %6.2f r(max)
assert r(max) < 4

* ---- numbers for the table: designs gamma = 0, 1, 2, 3, unknown ----
foreach g in 0 1 2 3 {
    foreach v in rej_lev_c rej_lev_w rej_lev_h3 rej_def_c rej_def_w releff {
        quietly summarize `v' if gamma == `g'
        local `v'_`g' : di %4.2f r(mean)
        local `v'_`g' = strtrim("``v'_`g''")
    }
}
foreach v in rej_lev_c rej_lev_w rej_lev_h3 rej_def_c rej_def_w releff {
    quietly summarize `v' if design == 10
    local `v'_u : di %4.2f r(mean)
    local `v'_u = strtrim("``v'_u'")
}

capture mkdir tables
file open tab using tables/tab_abk_heteroskedasticity_montecarlo.tex, write replace
file write tab "\begin{tabular}{@{}lccccccc@{}}" _n
file write tab "\toprule" _n
file write tab " & & \multicolumn{3}{c}{Levels OLS} & \multicolumn{2}{c}{Deflated by \$z\$} & \\" _n
file write tab "\cmidrule(lr){3-5}\cmidrule(lr){6-7}" _n
file write tab "Form of \$\sigma_i^2\$ & \$\gamma\$ & conv. & HC1 & HC3 & conv. & HC1 & s.d. ratio \\" _n
file write tab "\midrule" _n
file write tab "\$\sigma^2\$ & 0 & \textcolor{DeepNavy}{\bfseries `rej_lev_c_0'} & `rej_lev_w_0' & `rej_lev_h3_0' & \textcolor{Maroon}{\bfseries `rej_def_c_0'} & `rej_def_w_0' & `releff_0' \\" _n
file write tab "\$\sigma^2 z_i\$ & 1 & `rej_lev_c_1' & `rej_lev_w_1' & `rej_lev_h3_1' & \textcolor{Maroon}{\bfseries `rej_def_c_1'} & `rej_def_w_1' & `releff_1' \\" _n
file write tab "\$\sigma^2 z_i^2\$ & 2 & \textcolor{Maroon}{\bfseries `rej_lev_c_2'} & `rej_lev_w_2' & `rej_lev_h3_2' & \textcolor{DeepNavy}{\bfseries `rej_def_c_2'} & `rej_def_w_2' & \textcolor{DeepNavy}{\bfseries `releff_2'} \\" _n
file write tab "\$\sigma^2 z_i^3\$ & 3 & `rej_lev_c_3' & `rej_lev_w_3' & `rej_lev_h3_3' & \textcolor{Maroon}{\bfseries `rej_def_c_3'} & `rej_def_w_3' & `releff_3' \\" _n
file write tab "\$\sigma^2 z_i^2 (1 + 2D_i)^2\$ & unknown & `rej_lev_c_u' & `rej_lev_w_u' & `rej_lev_h3_u' & \textcolor{Maroon}{\bfseries `rej_def_c_u'} & \textcolor{DeepNavy}{\bfseries `rej_def_w_u'} & `releff_u' \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_abk_heteroskedasticity_montecarlo.tex
di "Done: tables/tab_abk_heteroskedasticity_montecarlo.tex"

* ---- figure 2: empirical size against gamma (log scale separates the lines) ----
keep if design <= 9
capture mkdir figures
twoway (line rej_lev_c gamma, lcolor(maroon) lwidth(thick)) ///
       (line rej_def_c gamma, lcolor(navy) lwidth(thick)) ///
       (line rej_lev_w gamma, lcolor(maroon) lwidth(medthick) lpattern(dash)) ///
       (line rej_def_w gamma, lcolor(navy) lwidth(medthick) lpattern(dash)), ///
    yline(0.05, lcolor(gs8) lpattern(dot)) ///
    xline(2, lcolor(gs10) lpattern(dash)) ///
    yscale(log) ///
    text(0.045 1.0 "nominal 5 percent (dotted)", color(gs6) size(small) placement(s)) ///
    text(`=rej_lev_c[8]' 3.5 "levels, conventional", color(maroon) size(medsmall) placement(n)) ///
    text(`=rej_def_c[8]' 3.5 "deflated, conventional", color(navy) size(medsmall) placement(s)) ///
    text(0.16 2.7 "levels, White HC1 (dashed)", color(maroon) size(medsmall) placement(n)) ///
    text(0.036 3.5 "deflated, White HC1 (dashed)", color(navy) size(medsmall) placement(s)) ///
    xtitle("Exponent {&gamma} in {&sigma}{sub:i}{sup:2} = {&sigma}{sup:2} z{sub:i}{sup:{&gamma}}") ///
    ytitle("Rejection rate of the true null (log scale)") ///
    xlabel(0(1)4) ylabel(0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1, angle(horizontal) format(%5.3f)) ///
    title("Each conventional t-test is correctly sized at one exponent") ///
    subtitle("1,000 samples of N = 740 per point (seed 1999); deflation by z assumes {&gamma} = 2") ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig2_abk_heteroskedasticity_size.png, replace width(2000)
di "Done: figures/fig2_abk_heteroskedasticity_size.png"

* ---- figure 3: relative efficiency of the deflated slope (log scale) ----
gen lreleff = log(releff)
twoway (line releff gamma, lcolor(navy) lwidth(thick)) ///
       (scatteri `=releff[5]' 2, mcolor(maroon) msymbol(circle) msize(large)), ///
    yline(1, lcolor(gs8) lpattern(dash)) ///
    xline(2, lcolor(gs10) lpattern(dash)) ///
    yscale(log) ///
    text(1.2 2.3 "deflated slope less precise than the levels slope", color(gs6) size(small) placement(e)) ///
    text(0.83 2.3 "deflated slope more precise than the levels slope", color(gs6) size(small) placement(e)) ///
    text(`=releff[5]' 2.1 "{&gamma} = 2: ratio `=strtrim(string(releff[5], "%4.2f"))'", color(maroon) size(medsmall) placement(e)) ///
    xtitle("Exponent {&gamma} in {&sigma}{sub:i}{sup:2} = {&sigma}{sup:2} z{sub:i}{sup:{&gamma}}") ///
    ytitle("s.d. of deflated slope over s.d. of levels slope (log scale)") ///
    xlabel(0(1)4) ylabel(0.02 0.05 0.1 0.2 0.5 1 2 5 10, angle(horizontal) format(%4.2f)) ///
    title("Weighting pays only when the weights match the form") ///
    subtitle("1,000 samples of N = 740 per point (seed 1999); a ratio below one favors deflation") ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig3_abk_heteroskedasticity_efficiency.png, replace width(2000)
di "Done: figures/fig3_abk_heteroskedasticity_efficiency.png"

log close
