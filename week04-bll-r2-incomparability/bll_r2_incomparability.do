********************************************************************************
* ACCTG 597A -- Week 4: R-Squared as a Sample Statistic
* Brown, Lo, and Lys (1999), "Use of R2 in Accounting Research: Measuring Changes
*   in Value Relevance over the Last Four Decades," Journal of Accounting and
*   Economics 28(2): 83-115.
* Gu (2007), "Across-Sample Incomparability of R2s and Additional Evidence on
*   Value Relevance Changes Over Time," Journal of Business Finance & Accounting
*   34(7-8): 1073-1098.
*
* Self-contained simulation, base Stata only (version 17 or later). Run from
* this folder:  do bll_r2_incomparability.do
* Figures are written to figures/ and booktabs table fragments to tables/.
* Runtime under a minute.
*
* Modules. In each, the true parameter is known by construction, so the
* property under study is visible directly.
*   1. Two samples, one relation (Gu 2007, pp. 1076-1077). x ~ N(2, 1) versus
*      x ~ N(2, 3), y = x + e with e ~ N(0, 1) in both. Same slope, same
*      residual s.d.; R2 = 0.50 versus 0.75. Identity: R2 = b^2 s2x/(b^2 s2x + s2e).
*   2. Three ways to move R2 (Gu's sample property; BLL's scale factor; a noisier
*      relation). Only the third changes the economic relation; the residual s.d.
*      separates it from the first but not from the second.
*   3. R2 against the scale factor's coefficient of variation (BLL Eq. (7)) at
*      alpha = 0 (monotone, Eq. (8')) and alpha = 12 (falls, then rises: the case
*      the 2002 erratum could not rule out). Identity: Eq. (7) equals
*      beta^2 V/(beta^2 V + sigma_e^2), V = sigma_w^2 + mu_w^2 omega_s^2/(1 + omega_s^2).
*   4. Gu (2007) Table 2 recreated: pricing errors rise less than in proportion
*      to scale, so a multiplicative (proportional) scale adjustment is
*      descriptively invalid.
********************************************************************************

version 17
clear all
set more off
capture log close
log using bll_r2_incomparability.log, replace text
capture mkdir figures
capture mkdir tables


********************************************************************************
* MODULE 1: two samples, one relation, two R2s (Gu 2007)
* DGP. Sample 1: x ~ N(2, 1). Sample 2: x ~ N(2, 3). Both: y = x + e, e ~ N(0, 1).
* N = 100,000 per sample. Population R2 = sigma_x^2/(sigma_x^2 + 1): 0.50 and 0.75.
* Identity checks (asserted to 1e-6): SST = b^2 Sxx + SSE; R2 = b^2 s2x/(b^2 s2x + s2e).
* Monte Carlo check (tolerance 0.02): R2 within 0.02 of the population value.
********************************************************************************
clear all
set seed 1999

local N    = 100000
local beta = 1
local sige = 1

forvalues k = 1/2 {
    clear
    set obs `N'
    if `k' == 1 local sdx = 1
    if `k' == 2 local sdx = sqrt(3)
    gen x = 2 + `sdx'*rnormal()
    gen e = `sige'*rnormal()
    gen y = `beta'*x + e
    quietly reg y x
    scalar b`k'   = _b[x]
    scalar se`k'  = _se[x]
    scalar r2_`k' = e(r2)
    scalar sse`k' = e(rss)
    scalar sst`k' = e(mss) + e(rss)
    scalar sde`k' = sqrt(e(rss)/(e(N) - 2))
    scalar s2e`k' = e(rss)/(e(N) - 1)
    quietly summarize x
    scalar sxx`k' = r(Var)*(r(N) - 1)
    scalar s2x`k' = r(Var)
    scalar sdx`k' = r(sd)
    scalar pop`k' = (`beta'^2*`sdx'^2)/(`beta'^2*`sdx'^2 + `sige'^2)
}

di _n "Module 1 identity checks (should be zero to machine precision):"
forvalues k = 1/2 {
    di "  Sample `k': [SST - (b^2 Sxx + SSE)]/SST       : " %12.8f ((sst`k' - (b`k'^2*sxx`k' + sse`k'))/sst`k')
    di "  Sample `k': R2 - b^2 s2x/(b^2 s2x + s2e)      : " %12.8f (r2_`k' - b`k'^2*s2x`k'/(b`k'^2*s2x`k' + s2e`k'))
    assert abs((sst`k' - (b`k'^2*sxx`k' + sse`k'))/sst`k') < 1e-6
    assert abs(r2_`k' - b`k'^2*s2x`k'/(b`k'^2*s2x`k' + s2e`k')) < 1e-6
}
di _n "Module 1 Monte Carlo checks (tolerance 0.02):"
forvalues k = 1/2 {
    di "  Sample `k': betahat = " %6.4f b`k' "  resid s.d. = " %6.4f sde`k' ///
       "  s.d.(x) = " %6.4f sdx`k' "  R2 = " %6.4f r2_`k' "  population R2 = " %6.4f pop`k'
    assert abs(r2_`k' - pop`k') < 0.02
}

forvalues k = 1/2 {
    local fb`k'   : di %4.2f b`k'
    local fb`k'   = strtrim("`fb`k''")
    local fse`k'  : di %5.3f se`k'
    local fse`k'  = strtrim("`fse`k''")
    local fsde`k' : di %4.2f sde`k'
    local fsde`k' = strtrim("`fsde`k''")
    local fsdx`k' : di %4.2f sdx`k'
    local fsdx`k' = strtrim("`fsdx`k''")
    local fr2`k'  : di %4.2f r2_`k'
    local fr2`k'  = strtrim("`fr2`k''")
    local fpop`k' : di %4.2f pop`k'
    local fpop`k' = strtrim("`fpop`k''")
}

file open tab using tables/tab_twosample.tex, write replace
file write tab "\begin{tabular}{@{}lcc@{}}" _n
file write tab "\toprule" _n
file write tab " & Sample 1 & Sample 2 \\" _n
file write tab "\midrule" _n
file write tab "\$\hat\beta\$ (s.e.) & `fb1' (`fse1') & `fb2' (`fse2') \\" _n
file write tab "\$\hat\sigma_\varepsilon\$ (residual s.d.) & \textcolor{DeepNavy}{\bfseries `fsde1'} & \textcolor{DeepNavy}{\bfseries `fsde2'} \\" _n
file write tab "\$\hat\sigma_x\$ & `fsdx1' & `fsdx2' \\" _n
file write tab "\midrule" _n
file write tab "\$R^2\$ & \textcolor{Maroon}{\bfseries `fr21'} & \textcolor{Maroon}{\bfseries `fr22'} \\" _n
file write tab "plim \$R^2\$ & `fpop1' & `fpop2' \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_twosample.tex
di "Done: tables/tab_twosample.tex"


********************************************************************************
* MODULE 2: three ways to move R2 (BLL Eqs. (1)-(2) with alpha = 0)
* DGP. w ~ N(2, sigma_w^2), eps ~ N(0, sigma_e^2), z = w + eps; s lognormal with
* E(s) = 1 and CV omega_s (s = 1 when omega_s = 0); y = s z, x = s w; regress y on x.
*   baseline: sigma_w = 1, omega_s = 0, sigma_e = 1
*   (a) sigma_w = sqrt(3);  (b) omega_s = 1;  (c) sigma_e = sqrt(3)
* Population values: R2 = V/(V + sigma_e^2), V = sigma_w^2 + 4 omega_s^2/(1 + omega_s^2);
*   resid s.d. = sqrt(1 + omega_s^2) sigma_e; s.d.(x) = sqrt(omega_s^2 sigma_w^2 + 4 omega_s^2 + sigma_w^2).
*   R2: 0.50, 0.75, 0.75, 0.25.  Resid s.d.: 1.00, 1.00, 1.41, 1.73.
* Identity check (asserted to 1e-6): R2 = b^2 s2x/(b^2 s2x + s2e) in each sample.
* Monte Carlo checks: R2 within 0.02, standard deviations within 0.05, of population values.
********************************************************************************
clear all
set seed 1999

local N    = 100000
local beta = 1
local muw  = 2

forvalues k = 0/3 {
    local sdw  = 1
    local ws   = 0
    local sige = 1
    if `k' == 1 local sdw  = sqrt(3)
    if `k' == 2 local ws   = 1
    if `k' == 3 local sige = sqrt(3)
    clear
    set obs `N'
    gen w   = `muw' + `sdw'*rnormal()
    gen eps = `sige'*rnormal()
    gen z   = `beta'*w + eps
    if `ws' == 0 {
        gen s = 1
    }
    else {
        local v = sqrt(ln(1 + `ws'^2))
        gen s = exp(-0.5*`v'^2 + `v'*rnormal())
    }
    gen y = s*z
    gen x = s*w
    quietly reg y x
    scalar b`k'   = _b[x]
    scalar r2_`k' = e(r2)
    scalar sde`k' = sqrt(e(rss)/(e(N) - 2))
    scalar s2e`k' = e(rss)/(e(N) - 1)
    quietly summarize x
    scalar s2x`k' = r(Var)
    scalar sdx`k' = r(sd)
    local V = `sdw'^2 + `muw'^2*`ws'^2/(1 + `ws'^2)
    scalar popr2`k'  = `beta'^2*`V'/(`beta'^2*`V' + `sige'^2)
    scalar popsde`k' = sqrt(1 + `ws'^2)*`sige'
    scalar popsdx`k' = sqrt(`ws'^2*`sdw'^2 + `ws'^2*`muw'^2 + `sdw'^2)
}

di _n "Module 2 identity checks (should be zero to machine precision):"
forvalues k = 0/3 {
    di "  Sample `k': R2 - b^2 s2x/(b^2 s2x + s2e) : " %12.8f (r2_`k' - b`k'^2*s2x`k'/(b`k'^2*s2x`k' + s2e`k'))
    assert abs(r2_`k' - b`k'^2*s2x`k'/(b`k'^2*s2x`k' + s2e`k')) < 1e-6
}
di _n "Module 2 Monte Carlo checks (R2 tolerance 0.02; s.d. tolerance 0.05):"
forvalues k = 0/3 {
    di "  Sample `k': R2 = " %6.4f r2_`k' " (pop " %6.4f popr2`k' ")  b = " %6.4f b`k' ///
       "  resid s.d. = " %6.4f sde`k' " (pop " %6.4f popsde`k' ")  s.d.(x) = " %6.4f sdx`k' " (pop " %6.4f popsdx`k' ")"
    assert abs(r2_`k' - popr2`k') < 0.02
    assert abs(sde`k' - popsde`k') < 0.05
    assert abs(sdx`k' - popsdx`k') < 0.05
}

forvalues k = 0/3 {
    local fb`k'   : di %4.2f b`k'
    local fb`k'   = strtrim("`fb`k''")
    local fr2`k'  : di %4.2f r2_`k'
    local fr2`k'  = strtrim("`fr2`k''")
    local fsde`k' : di %4.2f sde`k'
    local fsde`k' = strtrim("`fsde`k''")
    local fsdx`k' : di %4.2f sdx`k'
    local fsdx`k' = strtrim("`fsdx`k''")
}

file open tab using tables/tab_three_sources.tex, write replace
file write tab "\begin{tabular}{@{}lcccc@{}}" _n
file write tab "\toprule" _n
file write tab " & Baseline & (a) \$\sigma_w\$: 1 to 1.73 & (b) \$\omega_s\$: 0 to 1 & (c) \$\sigma_\varepsilon\$: 1 to 1.73 \\" _n
file write tab "\midrule" _n
file write tab "\$R^2\$ & `fr20' & \textcolor{Maroon}{\bfseries `fr21'} & \textcolor{Maroon}{\bfseries `fr22'} & `fr23' \\" _n
file write tab "\$\hat\beta\$ & `fb0' & `fb1' & `fb2' & `fb3' \\" _n
file write tab "\$\hat\sigma_\varepsilon\$ (residual s.d.) & `fsde0' & \textcolor{DeepNavy}{\bfseries `fsde1'} & `fsde2' & `fsde3' \\" _n
file write tab "\$\hat\sigma_x\$ & `fsdx0' & `fsdx1' & `fsdx2' & `fsdx3' \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_three_sources.tex
di "Done: tables/tab_three_sources.tex"


********************************************************************************
* MODULE 3: R2 against the CV of the scale factor (BLL Eq. (7); a figure)
* DGP. w ~ N(2, 1), eps ~ N(0, 1), beta = 1, z = alpha + w + eps; s lognormal with
* E(s) = 1 and CV omega_s on a grid from 0 to 2; N = 400,000 per grid point.
* Drawn curves: Eq. (7) on a fine grid (step 0.01). Markers: simulated R2.
* Identity check (asserted to 1e-10): Eq. (7) equals beta^2 V/(beta^2 V + sigma_e^2) at alpha = 0.
* Monte Carlo check (tolerance 0.02): simulated R2 within 0.02 of Eq. (7).
********************************************************************************
clear all
set seed 1999

local N    = 400000
local beta = 1
local muw  = 2
local sdw  = 1
local sige = 1

tempfile sw
tempname sims
postfile `sims' alpha ws r2sim r2eq7 r2v using `sw', replace
foreach a in 0 12 {
    foreach ws in 0 0.1 0.2 0.3 0.4 0.5 0.75 1 1.25 1.5 1.75 2 {
        clear
        set obs `N'
        gen w   = `muw' + `sdw'*rnormal()
        gen eps = `sige'*rnormal()
        gen z   = `a' + `beta'*w + eps
        if `ws' == 0 {
            gen s = 1
        }
        else {
            local v = sqrt(ln(1 + `ws'^2))
            gen s = exp(-0.5*`v'^2 + `v'*rnormal())
        }
        gen y = s*z
        gen x = s*w
        quietly reg y x
        local ws2 = `ws'^2
        local s2z = `beta'^2*`sdw'^2 + `sige'^2
        local muz = `a' + `beta'*`muw'
        local A   = `ws2'*`s2z' + `ws2'*`muz'^2 + `s2z'
        local D   = `ws2'*`sdw'^2 + `ws2'*`muw'^2 + `sdw'^2
        local c1  = `muw'*`ws2'/`D'
        local r2eq7 = `D'/`A'*(`beta' + `a'*`c1')^2
        local V   = `sdw'^2 + `muw'^2*`ws2'/(1 + `ws2')
        local r2v = `beta'^2*`V'/(`beta'^2*`V' + `sige'^2)
        post `sims' (`a') (`ws') (e(r2)) (`r2eq7') (`r2v')
    }
}
postclose `sims'
use `sw', clear
gen fine = 0

tempfile grid
preserve
clear
set obs 201
gen ws = (_n - 1)/100
expand 2
bysort ws: gen alpha = cond(_n == 1, 0, 12)
gen ws2 = ws^2
gen s2z = `beta'^2*`sdw'^2 + `sige'^2
gen muz = alpha + `beta'*`muw'
gen A   = ws2*s2z + ws2*muz^2 + s2z
gen D   = ws2*`sdw'^2 + ws2*`muw'^2 + `sdw'^2
gen c1  = `muw'*ws2/D
gen r2eq7 = D/A*(`beta' + alpha*c1)^2
gen fine = 1
keep alpha ws r2eq7 fine
save `grid'
restore
append using `grid'

gen gap_eq7_v = r2eq7 - r2v if alpha == 0 & fine == 0
gen gap_sim   = r2sim - r2eq7 if fine == 0
list alpha ws r2sim r2eq7 gap_sim if fine == 0, sep(12) noobs
di _n "Module 3 identity check, alpha = 0: max abs [Eq. (7) - beta^2 V/(beta^2 V + sigma_e^2)] (should be zero):"
quietly summarize gap_eq7_v
di "  " %14.10f max(abs(r(min)), abs(r(max)))
assert max(abs(r(min)), abs(r(max))) < 1e-10
di _n "Module 3 Monte Carlo check: max abs (simulated R2 - Eq. (7)) over 24 grid points (tolerance 0.02):"
quietly summarize gap_sim
di "  " %10.6f max(abs(r(min)), abs(r(max)))
assert max(abs(r(min)), abs(r(max))) < 0.02

quietly summarize r2eq7 if alpha == 0 & ws == 0 & fine == 1
local f00 : di %4.2f r(mean)
quietly summarize r2eq7 if alpha == 12 & fine == 1
local f8min : di %4.2f r(min)
local r8min = r(min)
quietly summarize ws if alpha == 12 & fine == 1 & abs(r2eq7 - `r8min') < 1e-9
local f8at : di %4.2f r(mean)
quietly summarize ws if alpha == 12 & fine == 1 & ws > 0 & r2eq7 >= 0.5
local f8cross : di %4.2f r(min)
local lim0 : di %4.2f (`beta'^2*(`sdw'^2 + `muw'^2))/(`beta'^2*`sdw'^2 + `sige'^2 + (`beta'*`muw')^2)
local Ewz8 = 12*`muw' + `beta'*(`sdw'^2 + `muw'^2)
local lim8 : di %4.2f (`Ewz8'^2)/((`sdw'^2 + `muw'^2)*(`beta'^2*`sdw'^2 + `sige'^2 + (12 + `beta'*`muw')^2))
di _n "Module 3 numbers: R2 at omega_s = 0: `f00'; alpha = 0 limit `lim0'; alpha = 12 minimum `f8min' at omega_s = `f8at', back to 0.50 at omega_s = `f8cross', limit `lim8'"

twoway (line r2eq7 ws if alpha == 0 & fine == 1, lcolor(navy) lwidth(thick))                      ///
       (scatter r2sim ws if alpha == 0 & fine == 0, mcolor(navy) msymbol(circle) msize(medium))   ///
       (line r2eq7 ws if alpha == 12 & fine == 1, lcolor(maroon) lwidth(thick))                    ///
       (scatter r2sim ws if alpha == 12 & fine == 0, mcolor(maroon) msymbol(circle_hollow) msize(medium)), ///
    yline(`f00', lcolor(gs8) lpattern(dash))                                           ///
    text(0.505 1.6 "scale-free R-squared = `f00'", color(gs6) size(medsmall) placement(n))        ///
    text(0.87 0.05 "{&alpha} = 0: rises toward `lim0'", color(navy) size(medium) placement(e))    ///
    text(0.37 0.55 "{&alpha} = 12: falls, then rises toward `lim8'", color(maroon) size(medium) placement(s)) ///
    xtitle("Coefficient of variation of the scale factor, {&omega}{sub:s}")            ///
    ytitle("R-squared of y = sz on x = sw")                                             ///
    xlabel(0(0.5)2) ylabel(0.3(0.1)0.9, angle(horizontal) format(%3.1f))               ///
    title("R-squared moves with the CV of scale, relation fixed")                       ///
    subtitle("Simulated data (seed 1999), N = 400,000 per point; lines are BLL Eq. (7)") ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig1_scale_sweep.png, replace width(2000)
di "Done: figures/fig1_scale_sweep.png"


********************************************************************************
* MODULE 4: Gu (2007) Table 2 recreated (levels model, pooled 1953-1998)
* The ten (mean |P-hat|, pricing error) pairs typed from Gu's Table 2, against a
* ray from the origin through decile 3 (what proportional scale would predict).
********************************************************************************
clear all

input decile scale perr
 1  3.444  5.285
 2  5.996  4.832
 3  7.870  5.613
 4  9.903  6.473
 5 12.178  7.440
 6 14.783  8.302
 7 18.024  8.968
 8 22.357  9.465
 9 28.795 10.468
10 45.369 13.439
end
scalar slope3      = 5.613/7.870
scalar ratio_scale = 45.369/7.870
scalar ratio_err   = 13.439/5.613
scalar ray10       = slope3*45.369
gen ray = slope3*scale
di _n "Module 4, recreated from Gu (2007), Table 2, levels model:"
di "  decile 10 to decile 3, scale ratio         : " %4.1f ratio_scale
di "  decile 10 to decile 3, pricing-error ratio : " %4.1f ratio_err
di "  proportional ray at decile 10 scale would predict a pricing error of " %5.1f ray10
local frs : di %3.1f ratio_scale
local fre : di %3.1f ratio_err

twoway (line ray scale, lcolor(maroon) lpattern(dash) lwidth(medthick))                  ///
       (connected perr scale, lcolor(navy) mcolor(navy) msymbol(circle) lwidth(thick) msize(medium)), ///
    text(30 30 "proportional to scale" "(ray through decile 3)", color(maroon) size(medsmall) placement(w) justification(right)) ///
    text(3.0 16 "Gu (2007) Table 2:" "pricing error by scale decile", color(navy) size(medsmall) placement(e) justification(left)) ///
    text(7.4 7.6 "decile 3", color(gs6) size(medsmall) placement(nw))                     ///
    text(8.3 45 "decile 10: scale `frs' to 1," "pricing error `fre' to 1", color(gs6) size(medsmall) placement(sw) justification(right)) ///
    xtitle("Mean absolute fitted price in the decile (dollars)")                          ///
    ytitle("Pricing error, root mean squared residual (dollars)")                        ///
    xlabel(0(10)50) ylabel(0(5)35, angle(horizontal))                                    ///
    title("Pricing errors rise less than in proportion to scale")                        ///
    subtitle("Recreated from Gu (2007) Table 2, levels model (2), pooled 1953-1998")     ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig2_gu_table2.png, replace width(2000)
di "Done: figures/fig2_gu_table2.png"

log close
