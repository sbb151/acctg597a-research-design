********************************************************************************
* ACCTG 597A -- Week 3 (Callen, Livnat, and Segal 2009 module)
* Controlling for a mediator: the credit rating between earnings and the
* CDS premium.  Pure simulation.
*
* Structural model (all shocks standard normal, scaled as shown):
*   rating   = -0.6 earn + 0.5 u + v          (ratings respond to earnings)
*   premium  = -0.4 earn + 1.0 rating + k u + w
* Total effect of earnings on the premium: -0.4 + (-0.6)(1.0) = -1.0.
* Direct effect (not through the rating): -0.4.
* Design A sets k = 0 (no unobserved credit shock shared by rating and
* premium); design B sets k = 0.5.  Each specification is estimated on
* 500 samples of N = 2,000; the figure plots the mean estimate.
* Without the rating control the regression recovers the total effect in
* both designs.  With the rating control it recovers the direct effect in
* design A and a biased quantity in design B, because conditioning on the
* rating induces a correlation between earnings and the shared shock.
* Requirements: base Stata. Output: figures/fig2_mediator.png
********************************************************************************

version 17
clear all
set more off
set seed 1363

local R = 500
local N = 2000

tempname mem
postfile `mem' rep bA_no bA_ctl bB_no bB_ctl using mediator_sims, replace

forvalues r = 1/`R' {
    quietly {
        clear
        set obs `N'
        gen earn = rnormal()
        gen u    = rnormal()
        gen v    = rnormal()
        gen w    = rnormal()
        gen rating = -0.6*earn + 0.5*u + v
        gen premA  = -0.4*earn + 1.0*rating + 0.0*u + w
        gen premB  = -0.4*earn + 1.0*rating + 0.5*u + w

        reg premA earn
        local bA_no = _b[earn]
        reg premA earn rating
        local bA_ctl = _b[earn]
        reg premB earn
        local bB_no = _b[earn]
        reg premB earn rating
        local bB_ctl = _b[earn]

        post `mem' (`r') (`bA_no') (`bA_ctl') (`bB_no') (`bB_ctl')
    }
}
postclose `mem'

use mediator_sims, clear
summarize bA_no bA_ctl bB_no bB_ctl
foreach v in bA_no bA_ctl bB_no bB_ctl {
    quietly summarize `v'
    local m_`v' = r(mean)
    local s_`v' = r(sd)
    di "`v': mean = " %7.3f `m_`v'' "   sd = " %6.3f `s_`v''
}

* analytical direct-effect estimand in design B (population regression of
* premB on earn and rating): the bias term equals 0.5 x cov(earn,u | rating)
* / var(earn | rating).
di "Design B: population coefficient on earn with rating control:"
di "  var(rating) = " 0.36 + 0.25 + 1
local vr = 0.36 + 0.25 + 1
local bias = 0.5 * (0.3/`vr') / (1 - 0.36/`vr')
di "  bias = " %6.3f `bias' "   coefficient = " %6.3f -0.4 + `bias'

clear
set obs 4
gen spec = _n
gen est = .
replace est = `m_bA_no'  in 1
replace est = `m_bA_ctl' in 2
replace est = `m_bB_no'  in 3
replace est = `m_bB_ctl' in 4
gen lab = string(est, "%5.2f")
gen ypos = est - 0.07

twoway (bar est spec if spec == 1 | spec == 3, barwidth(0.6) color(navy))    ///
       (bar est spec if spec == 2 | spec == 4, barwidth(0.6) color(maroon)),  ///
    yline(-1.0, lcolor(navy) lpattern(dash))                                  ///
    yline(-0.4, lcolor(maroon) lpattern(dash))                                ///
    yline(0, lcolor(gs6))                                                     ///
    text(-0.08 1 "`=string(`m_bA_no',"%5.2f")'", color(white) size(medsmall)) ///
    text(-0.08 2 "`=string(`m_bA_ctl',"%5.2f")'", color(white) size(medsmall)) ///
    text(-0.08 3 "`=string(`m_bB_no',"%5.2f")'", color(white) size(medsmall)) ///
    text(-0.08 4 "`=string(`m_bB_ctl',"%5.2f")'", color(white) size(medsmall)) ///
    text(-0.94 4.55 "Total effect = -1.0", color(navy) size(small) placement(w)) ///
    text(-0.46 4.55 "Direct effect = -0.4", color(maroon) size(small) placement(w)) ///
    xlabel(1 `""No rating" "control""' 2 `""Rating" "control""'              ///
           3 `""No rating" "control""' 4 `""Rating" "control""', noticks)     ///
    xscale(range(0.4 4.6))                                                    ///
    text(-1.33 1.5 "Design A: no shared shock", color(gs4) size(medsmall))    ///
    text(-1.33 3.5 "Design B: shared credit shock", color(gs4) size(medsmall)) ///
    ylabel(-1.2(0.2)0, angle(horizontal) format(%3.1f))                       ///
    yscale(range(-1.4 0.02))                                                  ///
    xtitle("") ytitle("Estimated coefficient on earnings", margin(r=2))                    ///
    title("Rating control removes the mediated effect") ///
    subtitle("Simulation; total effect -1.0, direct effect -0.4") ///
    legend(off)                                                               ///
    scheme(s1color)
graph export figures/fig2_mediator.png, replace width(2000)

capture erase mediator_sims.dta
di "Done: figures/fig2_mediator.png"
