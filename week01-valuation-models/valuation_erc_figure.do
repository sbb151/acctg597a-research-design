********************************************************************************
* ACCTG 597A -- Week 1 (Easton-Zmijewski module)
* Theoretical ERC as a function of earnings persistence.
*
* With AR(1) earnings persistence omega, discount rate r, and price equal to
* the present value of expected future earnings, the response of price to a
* $1 earnings innovation is
*     ERC(omega, r) = 1 + omega/(1 + r - omega),
* which nests the two certainty-world polar cases:
*     omega = 0 (purely transitory):  ERC = 1
*     omega = 1 (permanent):          ERC = (1+r)/r  (earnings capitalization)
* Requirements: base Stata. Output: figures/fig5_erc_persistence.png
********************************************************************************

version 17
clear all
set more off

clear
set obs 200
gen omega = (_n-1)/199
gen erc08 = 1 + omega/(1 + 0.08 - omega)
gen erc10 = 1 + omega/(1 + 0.10 - omega)
gen erc12 = 1 + omega/(1 + 0.12 - omega)

* shaded band: the range of typical estimated ERCs
gen lo = 1
gen hi = 3

twoway (rarea hi lo omega, color(gs14))                                       ///
       (line erc10 omega, lcolor(navy) lwidth(thick)),                        ///
    yline(11, lpattern(dot) lcolor(maroon))                                   ///
    text(12.2 0.08 "Permanent earnings (capitalization model): ERC = (1+r)/r = 11", color(maroon) size(small) place(e)) ///
    text(2.2 0.04 "Typical estimated ERCs: 1 to 3", color(gs6) size(small) place(e)) ///
    text(5.6 0.80 "ERC = 1 + {&omega}/(1+r{&minus}{&omega})", color(navy) size(small) place(w)) ///
    xtitle("Earnings persistence ({&omega})")                                 ///
    ytitle("Theoretical ERC (r = 10%)")                                       ///
    ylabel(0(2)14, angle(horizontal) format(%2.0f))                           ///
    xlabel(0(0.2)1, format(%2.1f))                                            ///
    title("The capitalization model predicts an ERC near 11")                 ///
    subtitle("Annual earnings are near a random walk, yet estimated ERCs are 1 to 3")  ///
    legend(off)                                                               ///
    scheme(s1color)
graph export figures/fig5_erc_persistence.png, replace width(2000)

di "Done: figures/fig5_erc_persistence.png"
