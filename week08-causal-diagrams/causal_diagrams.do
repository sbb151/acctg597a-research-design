********************************************************************************
* ACCTG 597A -- Week 8: Confounders, Mediators, Colliders, and the Arrows a Design Assumes Away
* Gow, I. D., D. F. Larcker, and P. C. Reiss (2016), "Causal Inference in
* Accounting Research," Journal of Accounting Research 54(2): 477-523.
*
* Self-contained simulation, base Stata only (version 17 or later). Run from
* this folder:  do causal_diagrams.do
* Figures are written to figures/ and booktabs table fragments to tables/.
* Runtime about one minute (the Monte Carlo in Module 3 dominates).
*
* Modules. In each, the true parameter is known by construction, so the
* property under study is visible directly. Throughout, the causal effect of
* x on y is beta = 0.30.
*   1. One third variable, three roles. The slope of y on x with and without z
*      as a control when z is a confounder (z -> x, z -> y), a mediator
*      (x -> z -> y), or a collider (x -> z <- y). Controlling for z removes the
*      bias in the first case, changes the estimand in the second, and creates a
*      bias in the third. Identity checks: the omitted-variable identity and the
*      mediation identity hold exactly; the collider slope matches its closed form.
*   2. Selection on a collider. Keeping only observations with a high value of
*      c = 0.7 x + 0.7 y + v (a variable caused by both x and y) biases the slope
*      of y on x in the retained sample, down to a negative slope in the top decile.
*   3. A small exclusion violation and a weak instrument. The instrument z has a
*      direct effect gamma = 0.05 on y (one sixth of beta). The population IV
*      slope is beta + gamma/pi, so as the first stage pi weakens the IV estimate
*      moves further from beta than OLS does, although OLS is confounded.
********************************************************************************

version 17
clear all
set more off
capture log close
log using causal_diagrams.log, replace text
capture mkdir figures
capture mkdir tables


********************************************************************************
* MODULE 1: one third variable, three roles
********************************************************************************
* DGP (N = 5,000, seed 2016); e and u are N(0, 0.8^2) unless stated.
*   Confounder: z ~ N(0,1); x = 0.6 z + u (so Var(x) = 1); y = 0.30 x + 0.50 z + e.
*               Slope of y on x = beta + gamma * delta, delta the slope of z on x
*               (0.6 in population): 0.30 + 0.50 * 0.6 = 0.60. Given z: 0.30.
*   Mediator:   x ~ N(0,1); m = 0.5 x + u; y = 0.10 x + 0.40 m + e.
*               Slope of y on x = total effect c' + a b = 0.10 + 0.5 * 0.4 = 0.30.
*               Given m: the direct effect c' = 0.10.
*   Collider:   x ~ N(0,1); y = 0.30 x + e; c = 0.7 x + 0.7 y + v, v ~ N(0, 0.5).
*               Slope of y on x = 0.30. Given c:
*               beta - dy (dx + dy beta) s2e / (dy^2 s2e + s2v) = -0.201.
clear all
set seed 2016
local N 5000
local beta 0.30

* ---- confounder ----
set obs `N'
gen double z = rnormal()
gen double x = 0.6*z + rnormal(0, 0.8)
gen double y = `beta'*x + 0.50*z + rnormal(0, 0.8)
quietly regress y x
scalar Zc_none  = _b[x]
scalar Zc_tnone = _b[x]/_se[x]
quietly regress y x z
scalar Zc_ctrl  = _b[x]
scalar Zc_tctrl = _b[x]/_se[x]
scalar Zc_gz    = _b[z]
quietly regress z x
scalar Zc_delta = _b[x]

* ---- mediator ----
clear
set obs `N'
gen double x = rnormal()
gen double m = 0.5*x + rnormal(0, 0.8)
gen double y = 0.10*x + 0.40*m + rnormal(0, 0.8)
quietly regress y x
scalar Zm_none  = _b[x]
scalar Zm_tnone = _b[x]/_se[x]
quietly regress y x m
scalar Zm_ctrl  = _b[x]
scalar Zm_tctrl = _b[x]/_se[x]
scalar Zm_b     = _b[m]
quietly regress m x
scalar Zm_a     = _b[x]

* ---- collider ----
clear
set obs `N'
gen double x = rnormal()
gen double e = rnormal(0, 0.8)
gen double y = `beta'*x + e
gen double c = 0.7*x + 0.7*y + rnormal(0, sqrt(0.5))
quietly regress y x
scalar Zk_none  = _b[x]
scalar Zk_tnone = _b[x]/_se[x]
quietly regress y x c
scalar Zk_ctrl  = _b[x]
scalar Zk_tctrl = _b[x]/_se[x]
scalar Zk_closed = `beta' - 0.7*(0.7 + 0.7*`beta')*0.64/(0.49*0.64 + 0.5)
tempfile collider
save `collider'

* ---- identity checks: print them for the log, then assert them ----
di _n "Module 1 identity checks (should be zero to machine precision):"
di "  confounder: b_short - (b_long + g_z * delta) : " %12.9f (Zc_none - (Zc_ctrl + Zc_gz*Zc_delta))
di "  mediator:   c - (c' + a b)                    : " %12.9f (Zm_none - (Zm_ctrl + Zm_a*Zm_b))
di "  collider:   slope given c minus closed form   : " %12.6f (Zk_ctrl - Zk_closed) "  (tolerance 0.03; sampling error only)"
assert abs(Zc_none - (Zc_ctrl + Zc_gz*Zc_delta)) < 1e-8
assert abs(Zm_none - (Zm_ctrl + Zm_a*Zm_b)) < 1e-8
assert abs(Zk_ctrl - Zk_closed) < 0.03

di _n "Module 1: slope of y on x (t-statistics in parentheses); the causal effect is 0.30"
di "  confounder: without z " %6.3f Zc_none " (" %4.1f Zc_tnone "), with z " %6.3f Zc_ctrl " (" %4.1f Zc_tctrl ")"
di "  mediator:   without m " %6.3f Zm_none " (" %4.1f Zm_tnone "), with m " %6.3f Zm_ctrl " (" %4.1f Zm_tctrl ")  [direct effect 0.10]"
di "  collider:   without c " %6.3f Zk_none " (" %4.1f Zk_tnone "), with c " %6.3f Zk_ctrl " (" %4.1f Zk_tctrl ")  [closed form " %6.3f Zk_closed "]"

* ---- numbers that travel to the table: formats that do not pad ----
foreach s in Zc_none Zc_ctrl Zm_none Zm_ctrl Zk_none Zk_ctrl {
    local f`s' : di %5.2f `s'
    local f`s' = strtrim("`f`s''")
}
foreach s in Zc_tnone Zc_tctrl Zm_tnone Zm_tctrl Zk_tnone Zk_tctrl {
    local f`s' : di %4.1f `s'
    local f`s' = strtrim("`f`s''")
}

* ---- booktabs fragment: escape math dollars as \$ (file write expands $) ----
file open tab using tables/tab_third_variable.tex, write replace
file write tab "\begin{tabular}{@{}llcc@{}}" _n
file write tab "\toprule" _n
file write tab "Role of \$z\$ & Arrows & \$y\$ on \$x\$ & \$y\$ on \$x\$ and \$z\$ \\" _n
file write tab "\midrule" _n
file write tab "Confounder & \$z \rightarrow x\$, \$z \rightarrow y\$ & \textcolor{Maroon}{\bfseries `fZc_none' (`fZc_tnone')} & \textcolor{DeepNavy}{\bfseries `fZc_ctrl' (`fZc_tctrl')} \\" _n
file write tab "Mediator & \$x \rightarrow z \rightarrow y\$ & \textcolor{DeepNavy}{\bfseries `fZm_none' (`fZm_tnone')} & `fZm_ctrl' (`fZm_tctrl') \\" _n
file write tab "Collider & \$x \rightarrow z \leftarrow y\$ & \textcolor{DeepNavy}{\bfseries `fZk_none' (`fZk_tnone')} & \textcolor{Maroon}{\bfseries `fZk_ctrl' (`fZk_tctrl')} \\" _n
file write tab "\midrule" _n
file write tab "\multicolumn{4}{@{}l}{Causal effect of \$x\$ on \$y\$: 0.30 in every design (direct effect 0.10 in the mediator design). Seed 2016, \$N\$ = 5,000.} \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_third_variable.tex
di "Done: tables/tab_third_variable.tex"


********************************************************************************
* MODULE 2: selection on a collider (a figure)
********************************************************************************
* The collider data of Module 1. Keep the observations with c above its
* (100 - q)th percentile, for q = 100, 90, ..., 10 percent retained, and
* estimate the slope of y on x in each retained sample.
use `collider', clear
tempname res
postfile `res' double(keep thr slope tstat n) using sel_results, replace
foreach q in 100 90 80 70 60 50 40 30 20 10 {
    if `q' == 100 {
        scalar Zthr = -99
    }
    else {
        quietly _pctile c, p(`=100-`q'')
        scalar Zthr = r(r1)
    }
    quietly regress y x if c >= Zthr
    post `res' (`q') (Zthr) (_b[x]) (_b[x]/_se[x]) (e(N))
}
postclose `res'
use sel_results, clear
di _n "Module 2: slope of y on x in the sample retained by c (the causal effect is 0.30)"
list keep slope tstat n, noobs sep(0)
local fs10 : di %5.2f slope[10]
local fs10 = strtrim("`fs10'")
twoway (line slope keep, lcolor(maroon) lwidth(thick))                                   ///
       (scatter slope keep, mcolor(maroon) msymbol(circle) msize(medium)),                ///
    yline(0.30, lcolor(navy) lpattern(dash) lwidth(medthick))                            ///
    yline(0, lcolor(gs10) lpattern(solid))                                                ///
    text(0.33 55 "causal effect 0.30", color(navy) size(medsmall) placement(n))          ///
    text(`=slope[10]-0.05' 14 "top decile: slope `fs10'", color(maroon) size(medsmall) placement(e)) ///
    xscale(reverse) xlabel(10(10)100) ylabel(-0.1(0.1)0.4, format(%4.1f))                 ///
    xtitle("Percent of the sample retained, by c (selection tightens to the right)")     ///
    ytitle("Slope of y on x in the retained sample")                                      ///
    title("Selecting on a collider bends the slope toward zero and below")               ///
    subtitle("Simulated data (seed 2016), N = 5,000; c = 0.7 x + 0.7 y + v")             ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig1_collider_selection.png, replace width(2000)
erase sel_results.dta
di "Done: figures/fig1_collider_selection.png"


********************************************************************************
* MODULE 3: a small exclusion violation and a weak instrument (a figure)
********************************************************************************
* DGP (seed 2016; N = 1,000 per draw; 1,000 draws per point):
*   z ~ N(0,1) the instrument; w ~ N(0,1) an unobserved confounder;
*   x = pi z + 0.5 w + v;  y = 0.30 x + 0.05 z + 0.5 w + e;  v, e ~ N(0,1).
* Population: OLS = 0.30 + (0.05 pi + 0.25)/(pi^2 + 1.25); IV = 0.30 + 0.05/pi.
* Identity (one draw, asserted): the 2SLS slope equals the reduced-form slope
* over the first-stage slope.
clear all
set seed 2016
local beta 0.30
local gamma 0.05
local N 1000
local R 1000

set obs `N'
gen double z = rnormal()
gen double w = rnormal()
gen double x = 0.3*z + 0.5*w + rnormal()
gen double y = `beta'*x + `gamma'*z + 0.5*w + rnormal()
quietly ivregress 2sls y (x = z)
scalar Zb_iv = _b[x]
quietly regress y z
scalar Zrf = _b[z]
quietly regress x z
scalar Zfs = _b[z]
di _n "Module 3 identity check (one draw, pi = 0.3; should be zero to machine precision):"
di "  b_2SLS - (reduced form / first stage) : " %12.9f (Zb_iv - Zrf/Zfs)
assert abs(Zb_iv - Zrf/Zfs) < 1e-8

mata:
rseed(2016)
beta = `beta'; gamma = `gamma'; N = `N'; R = `R'
pis = (1.0, 0.5, 0.3, 0.2, 0.1, 0.05)
out = J(0, 8, .)
for (k = 1; k <= cols(pis); k++) {
    pi = pis[k]
    biv = J(R, 1, .); bols = J(R, 1, .); F = J(R, 1, .)
    for (r = 1; r <= R; r++) {
        z = rnormal(N, 1, 0, 1); w = rnormal(N, 1, 0, 1)
        x = pi :* z :+ 0.5 :* w :+ rnormal(N, 1, 0, 1)
        y = beta :* x :+ gamma :* z :+ 0.5 :* w :+ rnormal(N, 1, 0, 1)
        X = (x, J(N, 1, 1)); Z = (z, J(N, 1, 1))
        bo = invsym(X'X) * X'y
        fsb = invsym(Z'Z) * Z'x
        rfb = invsym(Z'Z) * Z'y
        biv[r] = rfb[1] / fsb[1]
        bols[r] = bo[1]
        ex = x - Z*fsb
        s2 = (ex'ex) / (N - 2)
        F[r] = fsb[1]^2 / (s2 * invsym(Z'Z)[1,1])
    }
    s = sort(biv, 1)
    out = out \ (pi, beta + gamma/pi, s[floor(R/2)], s[floor(R/4)], s[floor(3*R/4)], mean(bols), beta + (gamma*pi + 0.25)/(pi^2 + 1.25), mean(F))
}
st_matrix("Zout", out)
end
matrix colnames Zout = pi popIV medIV q1IV q3IV meanOLS popOLS meanF
di _n "Module 3: 1,000 draws of N = 1,000 per point (the causal effect is 0.30)"
matrix list Zout, format(%8.3f)
assert abs(Zout[1,3] - (`beta' + `gamma')) < 0.02

clear
svmat double Zout, names(col)
gen double lpi = -log10(pi)
local fiv5 : di %4.2f medIV[5]
local fiv5 = strtrim("`fiv5'")
local fols5 : di %4.2f meanOLS[5]
local fols5 = strtrim("`fols5'")
twoway (rarea q1IV q3IV lpi, color(maroon%20) lwidth(none))                              ///
       (line popIV lpi, lcolor(maroon) lpattern(dash) lwidth(medium))                     ///
       (line medIV lpi, lcolor(maroon) lwidth(thick))                                     ///
       (scatter medIV lpi, mcolor(maroon) msymbol(circle))                                ///
       (line meanOLS lpi, lcolor(gs6) lpattern(shortdash) lwidth(medthick)),              ///
    yline(0.30, lcolor(navy) lpattern(dash) lwidth(medthick))                            ///
    text(0.26 0.05 "causal effect 0.30", color(navy) size(medsmall) placement(e))        ///
    text(`=meanOLS[1]+0.06' 0.05 "OLS, confounded by w", color(gs6) size(medsmall) placement(e)) ///
    text(`=medIV[5]+0.12' 0.75 "median IV `fiv5' at F = 9" "versus OLS `fols5'", color(maroon) size(medsmall) placement(w)) ///
    xlabel(0 "1 (F 804)" 0.301 "0.5 (F 204)" 0.523 "0.3 (F 74)" 0.699 "0.2 (F 33)" 1 "0.1 (F 9)" 1.301 "0.05 (F 3)", labsize(small)) ///
    ylabel(0(0.3)1.8, format(%3.1f))                                                      ///
    xtitle("First-stage coefficient (mean first-stage F); weaker instrument to the right") ///
    ytitle("Estimate of the effect of x on y")                                            ///
    title("A small exclusion violation grows as the instrument weakens")                   ///
    subtitle("1,000 draws of N = 1,000 (seed 2016); direct effect of z on y = 0.05; band = interquartile range") ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig2_iv_weak_instrument.png, replace width(2000)
di "Done: figures/fig2_iv_weak_instrument.png"

log close
