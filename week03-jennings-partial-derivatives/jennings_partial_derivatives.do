********************************************************************************
* ACCTG 597A -- Week 3: Coefficients as Partial Derivatives
* Jennings (1990, TAR), "A Note on Interpreting 'Incremental Information
* Content'," and the marginal-effects logic behind incremental-content tests.
*
* Self-contained simulation, base Stata only (version 17 or later). Run from
* this folder:  do jennings_partial_derivatives.do
* One figure is written to figures/ and three booktabs table fragments to
* tables/. Runtime a few seconds.
*
* Four modules. In each, a regression coefficient is a partial derivative
* taken while the other regressors are held fixed, and the marginal effect
* that answers the accounting question is a linear combination of
* coefficients (test it with lincom or test, not by reading one coefficient).
*   1. Jennings' identity. RET on (CFO, TA) versus RET on (CFO, INC) with
*      INC = CFO + TA: gamma1 = theta1 - theta2 and gamma2 = theta2 exactly.
*      When cash flow and accruals are valued equally, CFO given INC is zero
*      although CFO is strongly informative.
*   2. Interaction terms. In RET = b0 + b1 E + b2 Size + b3 E x Size + e,
*      b1 is the earnings response at Size = 0; the marginal effect is
*      b1 + b3 Size, and mean centering moves b1 to the sample mean.
*   3. Ohlson levels regression under clean surplus. In MV on E and BV_t with
*      BV_t = BV_t-1 + E - Div, a retained dollar of earnings moves value by
*      alpha1 + alpha2; regressing on E and K = BV_t-1 - Div returns it exactly.
*   4. BGL (1982) revisited. In the two-stage design the Panel (A) coefficient
*      on HC is the derivative along the first-stage line, b1 + c b2; holding
*      PRE fixed gives b1 and holding the depreciation adjustment fixed gives
*      b1 + b2. Same seed-228 data as week03-bgl-two-stage.
********************************************************************************

version 17
clear all
set more off
capture log close
log using jennings_partial_derivatives.log, replace text


********************************************************************************
* MODULE 1: Jennings identity, eq. (2) versus eq. (3)
********************************************************************************
clear all
set seed 1990

set obs 1000
gen cfo = 0.15*rnormal()
gen ta  = -0.4*cfo + 0.09*rnormal()
gen inc = cfo + ta
gen e   = 0.20*rnormal()
gen ret1 = 0.30*cfo + 0.30*ta + e
gen ret2 = 0.45*cfo + 0.15*ta + e

corr cfo ta

foreach k in 1 2 {
    * ---- eq. (2): components ----
    reg ret`k' cfo ta
    scalar th1_`k'  = _b[cfo]
    scalar th1t_`k' = _b[cfo]/_se[cfo]
    scalar th2_`k'  = _b[ta]
    scalar th2t_`k' = _b[ta]/_se[ta]
    scalar r2a_`k'  = e(r2)
    test cfo = ta
    scalar Feq_`k'  = r(F)
    * ---- eq. (3): component plus aggregate ----
    reg ret`k' cfo inc
    scalar g1_`k'   = _b[cfo]
    scalar g1t_`k'  = _b[cfo]/_se[cfo]
    scalar g2_`k'   = _b[inc]
    scalar g2t_`k'  = _b[inc]/_se[inc]
    scalar r2b_`k'  = e(r2)
    test cfo
    scalar Fzero_`k' = r(F)

    di _n "Economy `k' identity checks (should be zero to machine precision):"
    di "  gamma1 minus (theta1 - theta2) : " %10.6f (g1_`k' - (th1_`k' - th2_`k'))
    di "  gamma2 minus theta2            : " %10.6f (g2_`k' - th2_`k')
    di "  R2 eq.(2) minus R2 eq.(3)      : " %10.6f (r2a_`k' - r2b_`k')
    di "  F(theta1 = theta2) minus F(gamma1 = 0) : " %10.6f (Feq_`k' - Fzero_`k')
}

* ---- Write the booktabs fragment ----
foreach s in th1 th1t th2 th2t g1 g1t g2 g2t {
    foreach k in 1 2 {
        if substr("`s'", -1, 1) == "t" {
            local f`s'_`k' : di %4.1f `s'_`k'
        }
        else {
            local f`s'_`k' : di %5.2f `s'_`k'
        }
        local f`s'_`k' = strtrim("`f`s'_`k''")
    }
}
local fr2_1 : di %4.2f r2a_1
local fr2_1 = strtrim("`fr2_1'")
local fr2_2 : di %4.2f r2a_2
local fr2_2 = strtrim("`fr2_2'")

file open tab using tables/tab_reparam.tex, write replace
file write tab "\begin{tabular}{@{}lcccc@{}}" _n
file write tab "\toprule" _n
file write tab " & \multicolumn{2}{c}{Economy 1: \$\theta_1 = \theta_2 = 0.30\$} & \multicolumn{2}{c}{Economy 2: \$\theta_1 = 0.45,\ \theta_2 = 0.15\$} \\" _n
file write tab "\cmidrule(lr){2-3}\cmidrule(lr){4-5}" _n
file write tab " & Eq. (2) & Eq. (3) & Eq. (2) & Eq. (3) \\" _n
file write tab "\midrule" _n
file write tab "CFO & `fth1_1' (`fth1t_1') & \textcolor{Maroon}{\bfseries `fg1_1' (`fg1t_1')} & `fth1_2' (`fth1t_2') & \textcolor{Maroon}{\bfseries `fg1_2' (`fg1t_2')} \\" _n
file write tab "TA  & `fth2_1' (`fth2t_1') & & `fth2_2' (`fth2t_2') & \\" _n
file write tab "INC & & \textcolor{DeepNavy}{\bfseries `fg2_1' (`fg2t_1')} & & \textcolor{DeepNavy}{\bfseries `fg2_2' (`fg2t_2')} \\" _n
file write tab "\midrule" _n
file write tab "\$R^2\$ & `fr2_1' & `fr2_1' & `fr2_2' & `fr2_2' \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab

type tables/tab_reparam.tex
di "Done: tables/tab_reparam.tex"

********************************************************************************
* MODULE 2: interaction terms and mean centering
********************************************************************************
clear all
set seed 1990

set obs 1000
gen size = 6 + 1.5*rnormal()
gen earn = 0.08*rnormal()
gen ret  = 0.02 - 0.60*earn + 0.01*size + 0.35*earn*size + 0.20*rnormal()

* ---- Uncentered regression ----
gen earn_size = earn*size
reg ret earn size earn_size
scalar b1  = _b[earn]
scalar b1t = _b[earn]/_se[earn]
scalar b3  = _b[earn_size]
matrix V = e(V)
scalar v11 = V[1,1]
scalar v33 = V[3,3]
scalar v13 = V[1,3]

* ---- Mean-centered regression ----
sum size
scalar sbar = r(mean)
gen size_c = size - sbar
sum earn
gen earn_c = earn - r(mean)
gen earn_size_c = earn_c*size_c
reg ret earn_c size_c earn_size_c
scalar b1c  = _b[earn_c]
scalar b1ct = _b[earn_c]/_se[earn_c]
scalar b3c  = _b[earn_size_c]

di _n "Uncentered beta1 = " %6.3f b1 " (t = " %5.2f b1t "); beta3 = " %6.3f b3
di "Centered beta1   = " %6.3f b1c " (t = " %5.2f b1ct "); beta3 = " %6.3f b3c
di "Check: b1 + b3*mean(size) - b1c = " %10.6f (b1 + b3*sbar - b1c)
di "Check: b3 - b3c                  = " %10.6f (b3 - b3c)

* ---- Marginal effect line and 95 percent band over Size = 0 to 10 ----
_pctile size, p(2.5 97.5)
scalar slo = r(r1)
scalar shi = r(r2)

clear
set obs 101
gen s   = (_n - 1)/10
gen me  = b1 + b3*s
gen se  = sqrt(v11 + s^2*v33 + 2*s*v13)
gen lo  = me - 1.96*se
gen hi  = me + 1.96*se
gen inband = (s >= slo & s <= shi)

local fb1  : di %4.2f b1
local fb1t : di %4.1f b1t
local fb1c : di %4.2f b1c
local fb1ct: di %4.1f b1ct
local fsbar: di %3.1f sbar
local ys0 = b1
local ysc = b1c

twoway (rarea lo hi s if inband, color(navy%15) lwidth(none))                 ///
       (rarea lo hi s if !inband & s < sbar, color(gs14%40) lwidth(none))     ///
       (rarea lo hi s if !inband & s > sbar, color(gs14%40) lwidth(none))     ///
       (line me s, lcolor(navy) lwidth(thick))                                ///
       (scatteri `ys0' 0, mcolor(maroon) msymbol(circle) msize(large))       ///
       (scatteri `ysc' `=sbar', mcolor(navy) msymbol(circle) msize(large)),  ///
    yline(0, lcolor(gs8) lpattern(dash))                                     ///
    xline(`=slo', lcolor(gs10) lpattern(shortdash))                          ///
    xline(`=shi', lcolor(gs10) lpattern(shortdash))                          ///
    text(`=b1-0.6' 0.1 "uncentered {&beta}{sub:1} = `fb1'" "(t = `fb1t')" "the response at Size = 0", color(maroon) size(medium) placement(e) justification(left)) ///
    text(`=b1c+0.75' `=sbar-0.3' "centered {&beta}{sub:1} = `fb1c' (t = `fb1ct')" "the response at mean Size = `fsbar'", color(navy) size(medlarge) placement(w) justification(right)) ///
    text(-1.25 `=slo+0.15' "observed Size support (2.5th to 97.5th percentile)", color(gs6) size(medium) placement(e)) ///
    xtitle("Size (log market value)") ytitle("Earnings response: {&beta}{sub:1} + {&beta}{sub:3} Size") ///
    xlabel(0(2)10) ylabel(-1(1)3, angle(horizontal) format(%2.0f))            ///
    yscale(range(-1.4 3.2))                                                  ///
    title("The earnings response depends on where the derivative is evaluated") ///
    subtitle("Simulated data (seed 1990), N = 1,000; band is a 95 percent confidence interval") ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig1_interaction.png, replace width(2000)

di "Done: figures/fig1_interaction.png"

********************************************************************************
* MODULE 3: Ohlson levels regression under clean surplus
********************************************************************************
clear all
set seed 1995

set obs 1000
gen bv0 = 20 + 5*rnormal()
gen e   = 0.10*bv0 + 1.5*rnormal()
gen div = 0.4*e + 0.5*rnormal()
gen bv1 = bv0 + e - div
gen mv  = 2 + 5.0*e + 0.8*bv1 + 3*rnormal()

* ---- (1) MV on E and closing book value ----
reg mv e bv1
scalar a1   = _b[e]
scalar a1t  = _b[e]/_se[e]
scalar a2   = _b[bv1]
scalar a2t  = _b[bv1]/_se[bv1]
scalar r2_1 = e(r2)
lincom e + bv1
scalar tot  = r(estimate)
scalar tott = r(estimate)/r(se)

* ---- (2) MV on E and K = BV_t-1 - Div (= BV_t - E) ----
gen k = bv0 - div
reg mv e k
scalar c1   = _b[e]
scalar c1t  = _b[e]/_se[e]
scalar c2   = _b[k]
scalar c2t  = _b[k]/_se[k]
scalar r2_2 = e(r2)

di _n "Identity checks (should be zero to machine precision):"
di "  coef on E in (2) minus (alpha1 + alpha2) : " %10.6f (c1 - (a1 + a2))
di "  t on E in (2) minus lincom t             : " %10.6f (c1t - tott)
di "  coef on K in (2) minus alpha2            : " %10.6f (c2 - a2)
di "  R2 (1) minus R2 (2)                      : " %10.6f (r2_1 - r2_2)

* ---- Write the booktabs fragment ----
local fa1  : di %4.2f a1
local fa1 = strtrim("`fa1'")
local fa1t : di %4.1f a1t
local fa1t = strtrim("`fa1t'")
local fa2  : di %4.2f a2
local fa2 = strtrim("`fa2'")
local fa2t : di %4.1f a2t
local fa2t = strtrim("`fa2t'")
local fc1  : di %4.2f c1
local fc1 = strtrim("`fc1'")
local fc1t : di %4.1f c1t
local fc1t = strtrim("`fc1t'")
local fc2  : di %4.2f c2
local fc2 = strtrim("`fc2'")
local fc2t : di %4.1f c2t
local fc2t = strtrim("`fc2t'")
local fr1  : di %4.2f r2_1
local fr1 = strtrim("`fr1'")
local fr2  : di %4.2f r2_2
local fr2 = strtrim("`fr2'")

file open tab using tables/tab_ohlson.tex, write replace
file write tab "\begin{tabular}{@{}lcc@{}}" _n
file write tab "\toprule" _n
file write tab " & (1) & (2) \\" _n
file write tab " & \$MV\$ on \$E\$, \$BV_t\$ & \$MV\$ on \$E\$, \$K\$ \\" _n
file write tab "\midrule" _n
file write tab "\$E\$ & \textcolor{Maroon}{\bfseries `fa1' (`fa1t')} & \textcolor{DeepNavy}{\bfseries `fc1' (`fc1t')} \\" _n
file write tab "\$BV_t\$ & `fa2' (`fa2t') & \\" _n
file write tab "\$K = BV_{t-1} - Div\$ & & `fc2' (`fc2t') \\" _n
file write tab "\midrule" _n
file write tab "\$R^2\$ & `fr1' & `fr2' \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab

type tables/tab_ohlson.tex
di "Done: tables/tab_ohlson.tex"

********************************************************************************
* MODULE 4: BGL Panel (A) as a directional derivative
********************************************************************************
clear all
set seed 228

set obs 303
gen hc  = 0.134 + 0.216*rnormal()
gen pre = 0.148 + 1.46*(hc - 0.134) + 0.206*rnormal()
gen r   = 0.018 + 0.39*(hc - 0.134) - 0.02*(pre - 0.148) + 0.19*rnormal()

* ---- First stage (Panel (A)) ----
reg pre hc
scalar c = _b[hc]
predict zA, resid

* ---- One regression and the three directional derivatives ----
reg r hc pre
scalar b1 = _b[hc]
scalar b2 = _b[pre]
lincom hc
scalar d_pre  = r(estimate)
scalar d_pret = r(estimate)/r(se)
lincom hc + pre
scalar d_D    = r(estimate)
scalar d_Dt   = r(estimate)/r(se)
lincom hc + `=c'*pre
scalar d_Z    = r(estimate)
scalar d_Zt   = r(estimate)/r(se)

* ---- Cross-checks against the alternative parameterizations ----
reg r hc zA
scalar A_hc  = _b[hc]
scalar A_hct = _b[hc]/_se[hc]
gen d = pre - hc
reg r hc d
scalar I_hc  = _b[hc]
scalar I_hct = _b[hc]/_se[hc]
reg r hc
scalar biv_hc = _b[hc]

di _n "Identity checks (should be zero to machine precision):"
di "  lincom(hc + c*pre) minus Panel (A) HC coef : " %10.6f (d_Z - A_hc)
di "  its t minus Panel (A) HC t                 : " %10.6f (d_Zt - A_hct)
di "  Panel (A) HC coef minus bivariate slope    : " %10.6f (A_hc - biv_hc)
di "  lincom(hc + pre) minus coef on HC given D  : " %10.6f (d_D - I_hc)
di "  its t minus t on HC given D                : " %10.6f (d_Dt - I_hct)
di "  first-stage slope c = " %6.3f c

* ---- Write the booktabs fragment ----
local fpre  : di %4.2f d_pre
local fpre = strtrim("`fpre'")
local fpret : di %4.1f d_pret
local fpret = strtrim("`fpret'")
local fD    : di %4.2f d_D
local fD = strtrim("`fD'")
local fDt   : di %4.1f d_Dt
local fDt = strtrim("`fDt'")
local fZ    : di %4.2f d_Z
local fZ = strtrim("`fZ'")
local fZt   : di %4.1f d_Zt
local fZt = strtrim("`fZt'")
local fc    : di %4.2f c
local fc = strtrim("`fc'")

file open tab using tables/tab_bgl.tex, write replace
file write tab "\begin{tabular}{@{}llcc@{}}" _n
file write tab "\toprule" _n
file write tab "Held fixed & Derivative & Estimate (\$t\$) & Where it appears \\" _n
file write tab "\midrule" _n
file write tab "PRE & \$b_1\$ & \textcolor{Maroon}{\bfseries `fpre' (`fpret')} & Panel (B) residual; one regression \\" _n
file write tab "\$D = \mathrm{PRE} - \mathrm{HC}\$ & \$b_1 + b_2\$ & `fD' (`fDt') & HC in \$R\$ on HC, \$D\$ \\" _n
file write tab "\$Z = \mathrm{PRE} - a - c\,\mathrm{HC}\$ & \$b_1 + c\,b_2\$, \$c = `fc'\$ & \textcolor{DeepNavy}{\bfseries `fZ' (`fZt')} & Panel (A) HC; bivariate slope \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab

type tables/tab_bgl.tex
di "Done: tables/tab_bgl.tex"

log close
di "Done: all four modules."
