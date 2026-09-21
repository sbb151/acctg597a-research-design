********************************************************************************
* ACCTG 597A -- Week 5: The Mishkin Test as a Return Regression
* Sloan, R. G. 1996. Do stock prices fully reflect information in accruals and
* cash flows about future earnings? The Accounting Review 71 (3): 289-315.
*
* Self-contained simulation, base Stata only (version 17 or later). Run from
* this folder:  do sloan-mishkin-test.do
* Figures are written to figures/ and booktabs table fragments to tables/.
* Runtime: several minutes (Modules 2 and 4 are Monte Carlo experiments).
*
* The Mishkin (1983) test estimates a forecasting equation and a pricing
* equation jointly and asks whether the weights implied by prices (gamma*)
* equal the weights in the data (gamma). In each module the true parameters
* are known by construction, so the property under study is visible directly.
*   1. The unconstrained Mishkin estimates are an exact reparameterization of OLS.
*   2. The Mishkin LR statistic tracks the OLS Wald statistic (Monte Carlo).
*   3. An omitted, mispriced correlate makes rationally priced accruals appear mispriced.
*   4. Pooled LR statistics over-reject under cross-sectional correlation; the
*      Fama-MacBeth t-test does not.
*   5. Sloan's Tables 4 and 5 recreated as a figure.
*   6. The implied misperception equals -(return slope)/ERC (Lewellen 2010).
********************************************************************************

version 17
clear all
set more off
capture log close
log using sloan-mishkin-test.log, replace text
capture mkdir figures
capture mkdir tables

* ---------------------------------------------------------------------------
* program mishkin: the Mishkin (1983) test by stacked weighted nonlinear least
* squares.  Syntax: mishkin yfuture abret x1 x2 ... [, fast]
*   Forecasting eq.: yfuture = a0 + sum g_k x_k + v
*   Pricing eq.    : abret   = b*(yfuture - a0s - sum gs_k x_k) + e
* The two equations are stacked (2n rows) and each row is weighted by the
* inverse residual variance of its equation.  Option fast skips the nonlinear
* fit of the unconstrained system and uses its OLS equivalent (see Module 1).
* Returns r(LR), r(wald), r(b), r(g_k), r(gs_k) (Mishkin), r(go_k), r(gso_k) (OLS).
* ---------------------------------------------------------------------------
capture program drop mishkin
program define mishkin, rclass
    syntax varlist(min=3) [, fast]
    gettoken yf rest : varlist
    gettoken ar xs : rest
    local n = _N
    * OLS, equation by equation
    qui reg `yf' `xs'
    tempname s1 s2 bo
    scalar `s1' = e(rss)/e(N)
    local init "a0 `=_b[_cons]'"
    local k = 0
    foreach x of local xs {
        local ++k
        return scalar go_`k' = _b[`x']
        local init "`init' g`k' `=_b[`x']'"
    }
    qui reg `ar' `yf' `xs'
    scalar `s2' = e(rss)/e(N)
    scalar `bo' = _b[`yf']
    return scalar bo = `bo'
    local initu "`init' b `=`bo'' a0s `=-_b[_cons]/`bo''"
    local initc "`initu'"
    local k = 0
    foreach x of local xs {
        local ++k
        return scalar gso_`k' = -_b[`x']/`bo'
        local initu "`initu' gs`k' `=-_b[`x']/`bo''"
    }
    * OLS test that lagged information does not predict abnormal returns
    qui reg `ar' `xs'
    qui test `xs'
    return scalar wald = r(df)*r(F)
    * stack the two equations
    preserve
    tempvar id eq y d1 d2 w res
    gen long `id' = _n
    qui expand 2
    bysort `id': gen byte `eq' = _n
    gen double `y'  = cond(`eq'==1, `yf', `ar')
    gen byte `d1' = `eq'==1
    gen byte `d2' = `eq'==2
    gen double `w' = cond(`eq'==1, 1/`s1', 1/`s2')
    local fu "{a0}"
    local pu "{a0s}"
    local pc "{a0s}"
    local k = 0
    foreach x of local xs {
        local ++k
        local fu "`fu'+{g`k'}*`x'"
        local pu "`pu'+{gs`k'}*`x'"
        local pc "`pc'+{g`k'}*`x'"
    }
    if "`fast'"=="" {
        qui nl (`y' = `d1'*(`fu') + `d2'*{b}*(`yf'-(`pu'))) [aw=`w'], initial(`initu')
        return scalar b = _b[/b]
        local k = 0
        foreach x of local xs {
            local ++k
            return scalar g_`k'  = _b[/g`k']
            return scalar gs_`k' = _b[/gs`k']
        }
        qui predict double `res', resid
        qui replace `res' = `res'^2*`w'
        qui sum `res'
        local ssru = r(sum)
        drop `res'
    }
    else local ssru = 2*`n'
    qui nl (`y' = `d1'*(`fu') + `d2'*{b}*(`yf'-(`pc'))) [aw=`w'], initial(`initc')
    qui predict double `res', resid
    qui replace `res' = `res'^2*`w'
    qui sum `res'
    return scalar LR = 2*`n'*ln(r(sum)/`ssru')
    restore
end


********************************************************************************
* MODULE 1: The unconstrained Mishkin estimates are OLS
* DGP with Sloan (1996, Table 5) values as the truth; N = 40,679.
* The stacked weighted nonlinear system reproduces equation-by-equation OLS
* (gamma_k* = -delta_k/b, asserted to 1e-6), and LR is close to the OLS Wald statistic.
********************************************************************************
clear
set seed 1996
set obs 40679
gen acc = rnormal(0,0.08)
gen cf  = -0.5*acc + rnormal(0,0.08)
gen e1  = 0.01 + 0.765*acc + 0.855*cf + rnormal(0,0.06)
gen ar  = 1.894*(e1 - 0.01 - 0.911*acc - 0.826*cf) + rnormal(0,0.30)

mishkin e1 ar acc cf
di _n "Mishkin (stacked weighted NLS) versus OLS:"
di "  gamma1  : " %9.6f r(g_1)  "   OLS gamma1          : " %9.6f r(go_1)
di "  gamma2  : " %9.6f r(g_2)  "   OLS gamma2          : " %9.6f r(go_2)
di "  gamma1* : " %9.6f r(gs_1) "   OLS -delta1/b       : " %9.6f r(gso_1)
di "  gamma2* : " %9.6f r(gs_2) "   OLS -delta2/b       : " %9.6f r(gso_2)
di "  beta    : " %9.6f r(b)    "   OLS b               : " %9.6f r(bo)
di "  LR statistic: " %8.2f r(LR) "   OLS Wald chi2(2): " %8.2f r(wald)
di _n "Identity checks (should be zero to machine precision):"
di "  gamma1* minus (-delta1/b): " %12.8f (r(gs_1)-r(gso_1))
di "  gamma2* minus (-delta2/b): " %12.8f (r(gs_2)-r(gso_2))
assert abs(r(g_1)-r(go_1))   < 1e-6
assert abs(r(g_2)-r(go_2))   < 1e-6
assert abs(r(gs_1)-r(gso_1)) < 1e-6
assert abs(r(gs_2)-r(gso_2)) < 1e-6
assert abs(r(b)-r(bo))       < 1e-6
assert abs(r(LR)-r(wald))/r(wald) < 0.02

foreach s in g_1 g_2 gs_1 gs_2 b go_1 go_2 gso_1 gso_2 bo {
    local f`s' : di %5.3f r(`s')
    local f`s' = strtrim("`f`s''")
}
local fLR : di %6.1f r(LR)
local fLR = strtrim("`fLR'")
local fW  : di %6.1f r(wald)
local fW  = strtrim("`fW'")
file open tab using tables/tab_sloan_mishkin_test_identity.tex, write replace
file write tab "\begin{tabular}{@{}lccc@{}}" _n "\toprule" _n
file write tab " & Truth & Mishkin NLS & From OLS \\" _n "\midrule" _n
file write tab "\$\gamma_1\$ (accruals) & 0.765 & `fg_1' & `fgo_1' \\" _n
file write tab "\$\gamma_2\$ (cash flows) & 0.855 & `fg_2' & `fgo_2' \\" _n
file write tab "\$\gamma_1^{*}\$ & 0.911 & \textcolor{Maroon}{\bfseries `fgs_1'} & \textcolor{Maroon}{\bfseries `fgso_1'} \\" _n
file write tab "\$\gamma_2^{*}\$ & 0.826 & \textcolor{Maroon}{\bfseries `fgs_2'} & \textcolor{Maroon}{\bfseries `fgso_2'} \\" _n
file write tab "\$\beta\$ & 1.894 & `fb' & `fbo' \\" _n "\midrule" _n
file write tab "Test statistic, \$\chi^2(2)\$ & & LR \$=\$ `fLR' & Wald \$=\$ `fW' \\" _n
file write tab "\bottomrule" _n "\end{tabular}" _n
file close tab
type tables/tab_sloan_mishkin_test_identity.tex
di "Done: tables/tab_sloan_mishkin_test_identity.tex"


********************************************************************************
* MODULE 2: The LR statistic tracks the OLS Wald statistic
* Monte Carlo, 300 samples of N = 2,000, market accrual weight varied across samples.
********************************************************************************
clear
set seed 1996
tempname sims
tempfile mc
postfile `sims' gap LR wald using `mc', replace
forvalues r = 1/300 {
    qui {
        clear
        set obs 2000
        local gap = 0.30*runiform()
        gen acc = rnormal(0,0.08)
        gen cf  = -0.5*acc + rnormal(0,0.08)
        gen e1  = 0.01 + 0.765*acc + 0.855*cf + rnormal(0,0.06)
        gen ar  = 1.894*(e1 - 0.01 - (0.765+`gap')*acc - 0.855*cf) + rnormal(0,0.30)
        mishkin e1 ar acc cf, fast
        post `sims' (`gap') (r(LR)) (r(wald))
    }
}
postclose `sims'
use `mc', clear
corr LR wald
scalar rho = r(rho)
gen reldiff = abs(LR-wald)/wald
sum reldiff, meanonly
scalar mrd = r(mean)
gen byte dis = (LR>5.99) != (wald>5.99)
count if dis
scalar ndis = r(N)
di _n "Correlation of LR and Wald: " %8.5f rho
di "Mean relative difference  : " %8.5f mrd
di "Replications (of 300) in which the two tests disagree at 5 percent: " ndis
assert rho > 0.999
assert mrd < 0.02
local frho : di %7.5f rho
sum wald, meanonly
local mx = ceil(r(max)/10)*10
twoway (function y = x, range(0 `mx') lcolor(gs8) lpattern(dash))          ///
       (scatter LR wald, mcolor(navy%60) msymbol(circle) msize(small)),    ///
    xline(5.99, lcolor(maroon) lpattern(shortdash))                         ///
    yline(5.99, lcolor(maroon) lpattern(shortdash))                         ///
    text(`=0.55*`mx'' `=0.97*`mx'' "45-degree line", color(gs6) size(medium) placement(w)) ///
    text(`=0.90*`mx'' `=0.12*`mx'' "correlation = `frho'", color(navy) size(medium) placement(e)) ///
    text(`=0.12*`mx'' `=0.97*`mx'' "5 percent critical value of chi2(2), 5.99", color(maroon) size(medium) placement(w)) ///
    xtitle("OLS Wald statistic (returns on accruals and cash flows)")       ///
    ytitle("Mishkin LR statistic")                                          ///
    ylabel(, angle(horizontal))                                             ///
    title("The Mishkin LR statistic tracks the OLS Wald statistic")            ///
    subtitle("Simulated data (seed 1996): 300 samples, N = 2,000 each")     ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig2_sloan_mishkin_test_lrwald.png, replace width(2000)
di "Done: figures/fig2_sloan_mishkin_test_lrwald.png"


********************************************************************************
* MODULE 3: An omitted, mispriced correlate of accruals
* Accruals and cash flows are rationally priced; an omitted variable Z (correlated
* with accruals, mispriced) makes accruals appear mispriced. Analytical gap 0.146.
********************************************************************************
clear
set seed 1996
set obs 40679
gen acc = rnormal(0,0.08)
gen cf  = -0.5*acc + rnormal(0,0.08)
gen z   = 0.5*acc + rnormal(0,0.08)
gen e1  = 0.01 + 0.765*acc + 0.855*cf + 0*z + rnormal(0,0.06)
gen ar  = 1.894*(e1 - 0.01 - 0.765*acc - 0.855*cf - 0.292*z) + rnormal(0,0.30)

mishkin e1 ar acc cf
foreach s in g_1 gs_1 g_2 gs_2 {
    local o`s' : di %5.3f r(`s')
    local o`s' = strtrim("`o`s''")
}
local oLR : di %6.1f r(LR)
local oLR = strtrim("`oLR'")
scalar gap_o = r(gs_1)-r(g_1)
di _n "Z omitted : gamma1 = " %6.3f r(g_1) "  gamma1* = " %6.3f r(gs_1) "  gap = " %6.3f gap_o "  LR = " %7.1f r(LR)
di "            gamma2 = " %6.3f r(g_2) "  gamma2* = " %6.3f r(gs_2)
di "  analytical gap = 0.146 ; simulated minus analytical = " %7.4f (gap_o-0.146)
assert abs(gap_o-0.146) < 0.03

mishkin e1 ar acc cf z
foreach s in g_1 gs_1 g_2 gs_2 g_3 gs_3 {
    local i`s' : di %5.3f r(`s')
    local i`s' = strtrim("`i`s''")
}
local iLR : di %6.1f r(LR)
local iLR = strtrim("`iLR'")
scalar gap_i = r(gs_1)-r(g_1)
di "Z included: gamma1 = " %6.3f r(g_1) "  gamma1* = " %6.3f r(gs_1) "  gap = " %6.3f gap_i "  LR(3) = " %7.1f r(LR)
di "            gamma3 = " %6.3f r(g_3) "  gamma3* = " %6.3f r(gs_3)
assert abs(gap_i) < 0.03

file open tab using tables/tab_sloan_mishkin_test_omitted.tex, write replace
file write tab "\begin{tabular}{@{}lccc@{}}" _n "\toprule" _n
file write tab " & Truth & \$Z\$ omitted & \$Z\$ included \\" _n "\midrule" _n
file write tab "\$\gamma_1\$ (accruals, data) & 0.765 & `og_1' & `ig_1' \\" _n
file write tab "\$\gamma_1^{*}\$ (accruals, prices) & 0.765 & \textcolor{Maroon}{\bfseries `ogs_1'} & \textcolor{DeepNavy}{\bfseries `igs_1'} \\" _n
file write tab "\$\gamma_2\$ (cash flows, data) & 0.855 & `og_2' & `ig_2' \\" _n
file write tab "\$\gamma_2^{*}\$ (cash flows, prices) & 0.855 & `ogs_2' & `igs_2' \\" _n
file write tab "\$\gamma_3\$ (\$Z\$, data) & 0.000 & & `ig_3' \\" _n
file write tab "\$\gamma_3^{*}\$ (\$Z\$, prices) & 0.292 & & `igs_3' \\" _n "\midrule" _n
file write tab "LR statistic & & `oLR' {\small [\$\chi^2(2)\$]} & `iLR' {\small [\$\chi^2(3)\$]} \\" _n
file write tab "\bottomrule" _n "\end{tabular}" _n
file close tab
type tables/tab_sloan_mishkin_test_omitted.tex
di "Done: tables/tab_sloan_mishkin_test_omitted.tex"


********************************************************************************
* MODULE 4: Pooled LR statistics and cross-sectional correlation
* The null is true; the yearly return loading on accruals has mean zero. 500 replications.
* Runtime: a few minutes.
********************************************************************************
clear
set seed 1996
tempname sims
tempfile mc
postfile `sims' LR wald tfm using `mc', replace
forvalues r = 1/500 {
    qui {
        clear
        set obs 30
        gen year = _n
        gen ct = rnormal(0,0.40)
        expand 500
        gen acc = rnormal(0,0.08)
        gen cf  = -0.5*acc + rnormal(0,0.08)
        gen v   = rnormal(0,0.06)
        gen e1  = 0.01 + 0.765*acc + 0.855*cf + v
        gen ar  = 1.894*v + ct*acc + rnormal(0,0.30)
        mishkin e1 ar acc cf, fast
        local LR = r(LR)
        local W  = r(wald)
        gen bfm = .
        forvalues t = 1/30 {
            reg ar acc cf if year==`t'
            replace bfm = _b[acc] in `t'
        }
        sum bfm in 1/30
        local tfm = r(mean)/(r(sd)/sqrt(30))
        post `sims' (`LR') (`W') (`tfm')
    }
}
postclose `sims'
use `mc', clear
gen byte rLR = LR   > invchi2tail(2,0.05)
gen byte rW  = wald > invchi2tail(2,0.05)
gen byte rFM = abs(tfm) > invttail(29,0.025)
foreach s in rLR rW rFM {
    sum `s', meanonly
    scalar p_`s' = 100*r(mean)
    local f`s' : di %4.1f p_`s'
    local f`s' = strtrim("`f`s''")
}
di _n "Rejection rates at the nominal 5 percent level (null is true), 500 replications:"
di "  pooled Mishkin LR : " %5.1f p_rLR " percent"
di "  pooled OLS Wald   : " %5.1f p_rW  " percent"
di "  Fama-MacBeth t    : " %5.1f p_rFM " percent"
assert p_rLR > 20
assert abs(p_rFM - 5) < 3
file open tab using tables/tab_sloan_mishkin_test_size.tex, write replace
file write tab "\begin{tabular}{@{}lc@{}}" _n "\toprule" _n
file write tab "Test at the 5 percent level (a correct test rejects 5\%) & Rejection rate \\" _n "\midrule" _n
file write tab "Pooled Mishkin LR, \$\chi^2(2)\$ & \textcolor{Maroon}{\bfseries `frLR'\%} \\" _n
file write tab "Pooled OLS Wald, \$\chi^2(2)\$ & \textcolor{Maroon}{\bfseries `frW'\%} \\" _n
file write tab "Fama-MacBeth \$t\$ (one regression per year; \$t\$-test on 30 slopes) & \textcolor{DeepNavy}{\bfseries `frFM'\%} \\" _n
file write tab "\bottomrule" _n "\end{tabular}" _n
file close tab
type tables/tab_sloan_mishkin_test_size.tex
di "Done: tables/tab_sloan_mishkin_test_size.tex"


********************************************************************************
* MODULE 5: Sloan (1996), Tables 4 and 5, recreated as a figure
* No simulation.
********************************************************************************
clear
input x data price
1 0.841 0.840
2 0.765 0.911
3 0.855 0.826
end
gen xd = x - 0.08
gen xp = x + 0.08
twoway (pcspike data xd price xp, lcolor(gs10))                                   ///
       (scatter data xd, mcolor(navy) msymbol(circle) msize(vlarge))              ///
       (scatter price xp, mcolor(maroon) msymbol(diamond) msize(vlarge)),         ///
    yline(0.841, lcolor(gs8) lpattern(dash))                                       ///
    text(0.866 0.55 "fixation benchmark, 0.841", color(gs6) size(medsmall) placement(e)) ///
    text(0.765 1.86 "data: 0.765", color(navy) size(medium) placement(w))         ///
    text(0.911 2.14 "prices: 0.911", color(maroon) size(medium) placement(e))     ///
    text(0.870 2.86 "data: 0.855", color(navy) size(medium) placement(w))         ///
    text(0.812 3.14 "prices: 0.826", color(maroon) size(medium) placement(e))     ///
    text(0.815 1 "data 0.841, prices 0.840", color(gs4) size(medium))             ///
    text(0.740 1 "LR = 0.007 (p = 0.933)", color(gs4) size(medsmall))             ///
    text(0.740 2.5 "LR = 180.91 (p = 0.000)", color(gs4) size(medsmall))          ///
    xlabel(1 "Earnings" 2 "Accruals" 3 "Cash flows", noticks) xscale(range(0.5 3.95)) ///
    ylabel(0.70(0.05)0.95, angle(horizontal) format(%4.2f)) yscale(range(0.72 0.95)) ///
    xtitle("") ytitle("Persistence coefficient")                                  ///
    title("Prices overweight accruals and underweight cash flows")                ///
    subtitle("Recreated from Sloan (1996), Tables 4-5, pp. 304-305; N = 40,679")  ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig1_sloan_mishkin_test_table5.png, replace width(2000)
di "Done: figures/fig1_sloan_mishkin_test_table5.png"


********************************************************************************
* MODULE 6: The implied bias varies inversely with the ERC
* Recreated from Lewellen (2010, p. 461). No simulation.
********************************************************************************
clear
scalar gS = 1.894*(0.765-0.911)
scalar gD = 1.206*(0.647-0.938)
di "Return slope, Sloan: " %6.3f gS "   implied gap: " %6.3f (-gS/1.894)
di "Return slope, DRS  : " %6.3f gD "   implied gap: " %6.3f (-gD/1.206)
assert abs(-gS/1.894 - 0.146) < 1e-6
assert abs(-gD/1.206 - 0.291) < 1e-6
local gs = -gS
local gd = -gD
twoway (function y = `gs'/x, range(0.9 3) lcolor(navy) lwidth(medthick))            ///
       (function y = `gd'/x, range(0.9 3) lcolor(maroon) lwidth(medthick))          ///
       (scatteri 0.146 1.894, mcolor(navy) msymbol(circle) msize(vlarge))           ///
       (scatteri 0.291 1.206, mcolor(maroon) msymbol(diamond) msize(vlarge)),       ///
    text(0.085 1.87 "Sloan (1996):" "ERC 1.894, gap 0.146", color(navy) size(medium) placement(w) justification(right)) ///
    text(0.335 1.28 "Dechow, Richardson, and Sloan (2008):" "ERC 1.206, gap 0.291", color(maroon) size(medium) placement(e) justification(left)) ///
    text(0.062 2.95 "return slope -0.28", color(navy) size(medsmall) placement(w))  ///
    text(0.168 2.95 "return slope -0.35", color(maroon) size(medsmall) placement(w)) ///
    xtitle("Earnings response coefficient (beta)")                                   ///
    ytitle("Implied misperception (gamma1* minus gamma1)")                           ///
    xlabel(1(0.5)3) ylabel(0(0.1)0.4, angle(horizontal) format(%3.1f))               ///
    title("Similar return slopes imply different misperceptions")                   ///
    subtitle("Values from Lewellen (2010, p. 461); gap = -(return slope)/ERC")    ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig3_sloan_mishkin_test_erc.png, replace width(2000)
di "Done: figures/fig3_sloan_mishkin_test_erc.png"


log close
