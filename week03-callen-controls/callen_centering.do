********************************************************************************
* ACCTG 597A -- Week 3 (Callen, Livnat, and Segal 2009 module)
* Does demeaning a continuous regressor reduce the multicollinearity of its
* interaction term?  Pure simulation.
*
* Model: y = 1 + 0.5 x + 0.3 D + 0.4 (x*D) + e, with x ~ N(3, 1) so that x
* and x*D are highly collinear in raw form, D ~ Bernoulli(0.5), e ~ N(0, 1).
* For each of 500 samples of N = 500 we estimate the model with raw x and
* with demeaned x, and record (i) the t-statistic on the interaction and
* (ii) the variance inflation factor of the interaction regressor.
* Centering is a reparameterization: the column spaces coincide, so the
* interaction coefficient, its standard error, and its t-statistic are
* identical. The VIF falls, but only because the interaction regressor's
* own variance falls by the same factor.
* Requirements: base Stata. Output: figures/fig1_centering.png
********************************************************************************

version 17
clear all
set more off
set seed 597

local R = 500
local N = 500

tempname mem
postfile `mem' rep t_raw t_cen se_raw se_cen vif_raw vif_cen b_raw b_cen ///
    using centering_sims, replace

forvalues r = 1/`R' {
    quietly {
        clear
        set obs `N'
        gen x  = 3 + rnormal()
        gen D  = runiform() < 0.5
        gen xD = x*D
        gen y  = 1 + 0.5*x + 0.3*D + 0.4*xD + rnormal()

        * raw specification
        reg y x D xD
        local b_raw  = _b[xD]
        local se_raw = _se[xD]
        local t_raw  = _b[xD]/_se[xD]
        reg xD x D
        local vif_raw = 1/(1 - e(r2))

        * demeaned specification (x centered before forming the interaction)
        summarize x, meanonly
        gen xc  = x - r(mean)
        gen xcD = xc*D
        reg y xc D xcD
        local b_cen  = _b[xcD]
        local se_cen = _se[xcD]
        local t_cen  = _b[xcD]/_se[xcD]
        reg xcD xc D
        local vif_cen = 1/(1 - e(r2))

        post `mem' (`r') (`t_raw') (`t_cen') (`se_raw') (`se_cen') ///
            (`vif_raw') (`vif_cen') (`b_raw') (`b_cen')
    }
}
postclose `mem'

use centering_sims, clear
gen dt   = t_cen - t_raw
gen dse  = se_cen - se_raw
gen db   = b_cen - b_raw
summarize t_raw t_cen se_raw se_cen vif_raw vif_cen dt dse db
di "Max absolute difference in interaction t-statistics: " ///
    %12.10f max(abs(r(max)), abs(r(min)))
quietly summarize dt
local maxdt = max(abs(r(max)), abs(r(min)))
di "max |t_cen - t_raw| = " %12.10f `maxdt'
quietly summarize dse
di "max |se_cen - se_raw| = " %12.10f max(abs(r(max)), abs(r(min)))
quietly summarize vif_raw
local mvraw = r(mean)
quietly summarize vif_cen
local mvcen = r(mean)
di "Mean VIF of the interaction regressor: raw = " %5.2f `mvraw' ///
    ", centered = " %5.2f `mvcen'

* ---- left panel: interaction t-statistics, raw versus centered ----
twoway (function y = x, range(2 9) lcolor(gs8) lpattern(dash))               ///
       (scatter t_cen t_raw, mcolor(navy%60) msymbol(O) msize(small)),        ///
    xtitle("t-statistic on x*D, raw x")                                       ///
    ytitle("t-statistic on x*D, demeaned x", margin(r=2))                                  ///
    xlabel(2(2)8) ylabel(2(2)8, angle(horizontal))                            ///
    text(1.5 3.4 "On the 45-degree line",       ///
         color(gs6) size(small) placement(e))                                 ///
    title("t-statistic on x*D: unchanged", size(medium))                 ///
    legend(off) scheme(s1color) name(gt, replace)

* ---- right panel: VIF of the interaction regressor ----
twoway (function y = x, range(0 14) lcolor(gs8) lpattern(dash))              ///
       (scatter vif_cen vif_raw, mcolor(maroon%60) msymbol(D) msize(small)),  ///
    xtitle("VIF of x*D, raw x")                                               ///
    ytitle("VIF of x*D, demeaned x", margin(r=2))                                          ///
    xlabel(0(4)12) ylabel(0(4)12, angle(horizontal)) xscale(range(0 14))                          ///
    text(10.4 0.2 "VIF falls; SE identical",               ///
         color(gs6) size(small) placement(e))                                 ///
    title("VIF of x*D: falls", size(medium))                   ///
    legend(off) scheme(s1color) name(gv, replace)

graph combine gt gv, rows(1)                                                  ///
    title("Centering relabels; inference is unchanged") ///
    subtitle("Simulation: y = 1 + 0.5x + 0.3D + 0.4(x*D) + e, x ~ N(3,1), N = 500") ///
    scheme(s1color) iscale(1.15)
graph export figures/fig1_centering.png, replace width(2000)

capture erase centering_sims.dta
di "Done: figures/fig1_centering.png"
