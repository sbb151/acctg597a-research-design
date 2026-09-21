********************************************************************************
* ACCTG 597A -- Week 5: Parallel Trends as a Counterfactual Assumption
* Schafhautle and Veenman (2024), "Crowdsourced Forecasts and the Market
* Reaction to Earnings Announcement News," The Accounting Review 99(2): 421-456.
*
* Self-contained simulation, base Stata only (version 17 or later). Run from
* this folder:  do parallel_trends.do
* Figures are written to figures/ and booktabs table fragments to tables/.
* Runtime about five minutes (Module 2 is a Monte Carlo experiment with
* 14,000 simulated panels).
*
* Modules. In each, the true treatment effect is known by construction, and
* the untreated potential outcome Y(0) of treated firms after treatment, which
* no real dataset records, is available to the simulator.
*   1. Three worlds. The pre-trend test passes although parallel trends fails
*      (World 2, whose observed panel is identical to World 1) and rejects
*      although parallel trends holds (World 3, a transitory pre-period shock).
*   2. The pre-trend test as a screen. With a differential drift in Y(0) and a
*      true effect of zero, the pre-trend test rejects far less often than the
*      DiD t-test rejects the true null.
********************************************************************************

version 17
clear all
set more off
capture log close
log using parallel_trends.log, replace text
capture mkdir figures
capture mkdir tables


********************************************************************************
* MODULE 1: three worlds
*
*
* What it demonstrates. The PTA restricts the mean untreated potential outcome
* Y(0) of the treated group after treatment, which is never observed. A
* pre-trend test is computed from observed outcomes before treatment.
*   World 1: the PTA holds and the pre-trend test passes (the benchmark).
*   World 2: the PTA FAILS and the pre-trend test PASSES. The observed panel is
*            identical to World 1 row by row, so no statistic can tell them
*            apart.
*   World 3: the PTA HOLDS and the pre-trend test REJECTS. Treated firms
*            experience a transitory shock at k = -3 that reverses at k = -2
*            (e.g., an accrual that reverses). The pre-period gaps average
*            zero, so both the pooled coefficient on Post and the event-study
*            coefficients relative to k = -1 are unbiased for the ATT.
*
* DGP. N = 4,000 firms (2,000 treated), event years k = -3,...,+3, treatment at
* k = 0 for all treated firms. alpha_i ~ N(D_i, 1); lambda_k = 0.1*k;
* eps_ik ~ N(0, 1). path_k = (0.2, 0.4, 0.5, 0.5) for k = 0, 1, 2, 3.
* blip_k = (-0.25, +0.25) for k = -3, -2 and zero from k = -1 onward.
*   World 1: Y(0) = alpha + lambda + eps;               Y(1) = Y(0) + path_k.
*   World 2: Y(0) = alpha + lambda + eps + D*path_k;    Y(1) = Y(0).
*   World 3: Y(0) = alpha + lambda + eps + D*blip_k;    Y(1) = Y(0) + path_k.
*   Observed: Y = Y(1) if D = 1 and k >= 0, and Y = Y(0) otherwise.
*
* Identity checks (machine precision):
*   (a) observed Y is the same in Worlds 1 and 2, row by row;
*   (b) each event-study coefficient equals the double difference of the four
*       group-by-period means relative to k = -1;
*   (c) beta_k(Y) = beta_k(Y(0)) + beta_k(Y - Y(0)), in Worlds 1 and 2;
*   (d) with a balanced panel and one treatment date, the two-way fixed
*       effects coefficient on Post equals the 2x2 DiD of means;
*   (e) beta_k(World 3) - beta_k(World 1) = blip_k: the worlds differ before
*       treatment only, and the post-period coefficients coincide;
*   (f) the pooled coefficient on Post is the same in Worlds 1 and 3, because
*       the pre-period gaps of World 3 average zero.
*
* Estimation: xtreg, fe with firm-clustered standard errors (areg would count
* the absorbed firm effects in the cluster correction and overstate them).
*
* Requirements: base Stata. Output: tables/tab_sv_parallel_trends_worlds.tex,
* figures/fig1a_sv_parallel_trends_w12.png, figures/fig1b_sv_parallel_trends_w3.png
********************************************************************************
clear all
set seed 2024

* ---- the panel ----
set obs 4000
gen firm = _n
gen D = firm <= 2000
gen double alpha = D + rnormal()
expand 7
bysort firm: gen t = _n
gen k = t - 4
gen post = k >= 0
gen double eps = rnormal()
gen double path  = 0.2*(k==0) + 0.4*(k==1) + 0.5*(k==2) + 0.5*(k==3)
gen double blip  = -0.25*(k==-3) + 0.25*(k==-2)

* ---- potential outcomes and observed outcomes, world by world ----
gen double y0_w1 = alpha + 0.1*k + eps
gen double y1_w1 = y0_w1 + path
gen double y_w1  = cond(D==1 & post==1, y1_w1, y0_w1)

gen double y0_w2 = y0_w1 + D*path
gen double y1_w2 = y0_w2
gen double y_w2  = cond(D==1 & post==1, y1_w2, y0_w2)

gen double y0_w3 = y0_w1 + D*blip
gen double y1_w3 = y0_w3 + path
gen double y_w3  = cond(D==1 & post==1, y1_w3, y0_w3)

* ---- identity (a): Worlds 1 and 2 share one observed dataset ----
gen double gap = abs(y_w1 - y_w2)
quietly summarize gap
di _n "Identity checks (should be zero to machine precision):"
di "  (a) max abs difference in observed Y, Worlds 1 and 2   : " %12.8f r(max)
assert r(max) < 1e-6

* ---- event-time indicators, k = -1 omitted ----
forvalues j = 1/7 {
    gen dk`j' = D*(t==`j')
}
drop dk3
gen DPost = D*post
xtset firm t

* ---- event studies: observed Y, and the infeasible regressions on Y(0) ----
foreach v in y_w1 y_w3 y0_w1 y0_w2 y0_w3 {
    quietly xtreg `v' i.t dk1 dk2 dk4 dk5 dk6 dk7, fe vce(cluster firm)
    foreach j in 1 2 4 5 6 7 {
        scalar b_`v'_`j'  = _b[dk`j']
        scalar se_`v'_`j' = _se[dk`j']
    }
    quietly testparm dk1 dk2
    scalar pF_`v' = r(p)
    scalar F_`v'  = r(F)
    quietly lincom (dk4 + dk5 + dk6 + dk7)/4
    scalar avg_`v'  = r(estimate)
    scalar avgt_`v' = r(estimate)/r(se)
}
* the treatment-effect component Y - Y(0), world 1 (it is zero in world 2)
gen double te_w1 = y_w1 - y0_w1
quietly xtreg te_w1 i.t dk1 dk2 dk4 dk5 dk6 dk7, fe
foreach j in 1 2 4 5 6 7 {
    scalar b_te_`j' = _b[dk`j']
}

* ---- identity (b): coefficients are double differences of group means ----
scalar maxb = 0
foreach j in 1 2 4 5 6 7 {
    quietly summarize y_w1 if D==1 & t==`j'
    scalar m1j = r(mean)
    quietly summarize y_w1 if D==1 & t==3
    scalar m1b = r(mean)
    quietly summarize y_w1 if D==0 & t==`j'
    scalar m0j = r(mean)
    quietly summarize y_w1 if D==0 & t==3
    scalar m0b = r(mean)
    scalar maxb = max(maxb, abs(b_y_w1_`j' - ((m1j - m1b) - (m0j - m0b))))
}
di "  (b) max abs(beta_k minus double difference of means)   : " %12.8f maxb
assert maxb < 1e-6

* ---- identity (c): beta_k(Y) = beta_k(Y(0)) + beta_k(Y - Y(0)) ----
scalar maxc1 = 0
scalar maxc2 = 0
foreach j in 1 2 4 5 6 7 {
    scalar maxc1 = max(maxc1, abs(b_y_w1_`j' - (b_y0_w1_`j' + b_te_`j')))
    scalar maxc2 = max(maxc2, abs(b_y_w1_`j' - b_y0_w2_`j'))
}
di "  (c) world 1: max abs(beta(Y) - beta(Y0) - beta(Y-Y0))  : " %12.8f maxc1
di "  (c) world 2: max abs(beta(Y) - beta(Y0))               : " %12.8f maxc2
assert maxc1 < 1e-6
assert maxc2 < 1e-6

* ---- pooled two-way fixed effects (eq. 1) and identity (d) ----
foreach v in y_w1 y_w3 y0_w1 y0_w2 y0_w3 {
    quietly xtreg `v' i.t DPost, fe vce(cluster firm)
    scalar d_`v'  = _b[DPost]
    scalar dt_`v' = _b[DPost]/_se[DPost]
}
quietly summarize y_w1 if D==1 & post==1
scalar a11 = r(mean)
quietly summarize y_w1 if D==1 & post==0
scalar a10 = r(mean)
quietly summarize y_w1 if D==0 & post==1
scalar a01 = r(mean)
quietly summarize y_w1 if D==0 & post==0
scalar a00 = r(mean)
di "  (d) TWFE delta minus 2x2 DiD of means                  : " %12.8f (d_y_w1 - ((a11-a10)-(a01-a00)))
assert abs(d_y_w1 - ((a11-a10)-(a01-a00))) < 1e-6

* ---- identity (e): World 3 differs from World 1 before treatment only ----
scalar maxe = 0
scalar maxe = max(maxe, abs((b_y_w3_1 - b_y_w1_1) - (-0.25)))
scalar maxe = max(maxe, abs((b_y_w3_2 - b_y_w1_2) - (0.25)))
foreach j in 4 5 6 7 {
    scalar maxe = max(maxe, abs(b_y_w3_`j' - b_y_w1_`j'))
}
di "  (e) max abs(beta_k(W3) - beta_k(W1) - blip_k)          : " %12.8f maxe
assert maxe < 1e-6
di "  (f) pooled delta, World 3 minus World 1                : " %12.8f (d_y_w3 - d_y_w1)
assert abs(d_y_w3 - d_y_w1) < 1e-6

di _n "Numbers for the slides:"
di "  pooled delta on Post, observed Y, Worlds 1 and 2   : " %7.4f d_y_w1 "  (t = " %5.2f dt_y_w1 ")"
di "  pooled delta on Post, observed Y, World 3          : " %7.4f d_y_w3 "  (t = " %5.2f dt_y_w3 ")"
di "  pre-trend F test p-value, Worlds 1 and 2           : " %7.4f pF_y_w1 "  (F = " %6.2f F_y_w1 ")"
di "  pre-trend F test p-value, World 3                  : " %7.4f pF_y_w3 "  (F = " %6.2f F_y_w3 ")"
di "  infeasible pooled delta on Y(0), World 1           : " %7.4f d_y0_w1 "  (t = " %5.2f dt_y0_w1 ")"
di "  infeasible pooled delta on Y(0), World 2           : " %7.4f d_y0_w2 "  (t = " %5.2f dt_y0_w2 ")"
di "  infeasible pooled delta on Y(0), World 3           : " %7.4f d_y0_w3 "  (t = " %5.2f dt_y0_w3 ")"
di "  mean of beta_0..beta_3, observed Y, all worlds     : " %7.4f avg_y_w1 "  (t = " %5.2f avgt_y_w1 ")"
di "  true average ATT, Worlds 1 / 2 / 3                 :  0.4000 / 0.0000 / 0.4000"
foreach j in 1 2 4 5 6 7 {
    di "  k = " %2.0f (`j'-4) ": b(Y,w1=w2) = " %6.3f b_y_w1_`j' "  b(Y,w3) = " %6.3f b_y_w3_`j' "  b(Y0,w1) = " %6.3f b_y0_w1_`j' "  b(Y0,w2) = " %6.3f b_y0_w2_`j' "  b(Y0,w3) = " %6.3f b_y0_w3_`j'
}

* ---- numbers that travel to the slide ----
foreach s in d_y_w1 d_y_w3 pF_y_w1 {
    local f`s' : di %4.2f `s'
    local f`s' = strtrim("`f`s''")
}
foreach s in dt_y_w1 dt_y_w3 {
    local f`s' : di %3.1f `s'
    local f`s' = strtrim("`f`s''")
}
local bias1 : di %5.2f (d_y_w1 - 0.40)
local bias1 = strtrim("`bias1'")
local bias2 : di %5.2f (d_y_w1 - 0.00)
local bias2 = strtrim("`bias2'")
local bias3 : di %5.2f (d_y_w3 - 0.40)
local bias3 = strtrim("`bias3'")
assert pF_y_w3 < 0.001

* ---- booktabs fragment ----
file open tab using tables/tab_sv_parallel_trends_worlds.tex, write replace
file write tab "\begin{tabular}{@{}lccc@{}}" _n
file write tab "\toprule" _n
file write tab " & World 1 & World 2 & World 3 \\" _n
file write tab "\midrule" _n
file write tab "\multicolumn{4}{@{}l}{\emph{Known only to the simulator}} \\" _n
file write tab "\quad Parallel trends, \$k \geq 0\$ & holds & \textcolor{Maroon}{\bfseries fails} & holds \\" _n
file write tab "\quad Average ATT & 0.40 & \textcolor{Maroon}{\bfseries 0.00} & 0.40 \\" _n
file write tab "\midrule" _n
file write tab "\multicolumn{4}{@{}l}{\emph{Computed from observed outcomes}} \\" _n
file write tab "\quad Pre-trend \$F\$ test, \$p\$ & `fpF_y_w1' & \textcolor{Maroon}{\bfseries `fpF_y_w1'} & \textcolor{Maroon}{\bfseries \$<\$\,0.001} \\" _n
file write tab "\quad \$\hat\delta\$ on \$Post\$ (\$t\$) & `fd_y_w1' (`fdt_y_w1') & `fd_y_w1' (`fdt_y_w1') & `fd_y_w3' (`fdt_y_w3') \\" _n
file write tab "\quad \$\hat\delta\$ minus ATT & \$`bias1'\$ & \textcolor{Maroon}{\bfseries `bias2'} & \$`bias3'\$ \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_sv_parallel_trends_worlds.tex

* ---- figures: coefficients into a small dataset ----
preserve
clear
set obs 7
gen k = _n - 4
foreach s in obs12 obs3 w1 w2 w3 {
    gen b_`s' = 0
    gen lo_`s' = .
    gen hi_`s' = .
}
foreach j in 1 2 4 5 6 7 {
    replace b_obs12  = b_y_w1_`j'                       in `j'
    replace lo_obs12 = b_y_w1_`j' - 1.96*se_y_w1_`j'    in `j'
    replace hi_obs12 = b_y_w1_`j' + 1.96*se_y_w1_`j'    in `j'
    replace b_obs3   = b_y_w3_`j'                       in `j'
    replace lo_obs3  = b_y_w3_`j' - 1.96*se_y_w3_`j'    in `j'
    replace hi_obs3  = b_y_w3_`j' + 1.96*se_y_w3_`j'    in `j'
    replace b_w1     = b_y0_w1_`j'                      in `j'
    replace lo_w1    = b_y0_w1_`j' - 1.96*se_y0_w1_`j'  in `j'
    replace hi_w1    = b_y0_w1_`j' + 1.96*se_y0_w1_`j'  in `j'
    replace b_w2     = b_y0_w2_`j'                      in `j'
    replace lo_w2    = b_y0_w2_`j' - 1.96*se_y0_w2_`j'  in `j'
    replace hi_w2    = b_y0_w2_`j' + 1.96*se_y0_w2_`j'  in `j'
    replace b_w3     = b_y0_w3_`j'                      in `j'
    replace lo_w3    = b_y0_w3_`j' - 1.96*se_y0_w3_`j'  in `j'
    replace hi_w3    = b_y0_w3_`j' + 1.96*se_y0_w3_`j'  in `j'
}
gen k_l = k - 0.22
gen k_r = k + 0.22

* Figure 1a: Worlds 1 and 2 (the test passes; the PTA holds in one, fails in the other)
twoway (rcap lo_w1 hi_w1 k_l, lcolor(gs8))                                       ///
       (scatter b_w1 k_l, mcolor(gs6) msymbol(circle_hollow) msize(large))       ///
       (rcap lo_w2 hi_w2 k_r, lcolor(maroon))                                    ///
       (scatter b_w2 k_r, mcolor(maroon) msymbol(diamond) msize(large))          ///
       (rcap lo_obs12 hi_obs12 k, lcolor(navy))                                  ///
       (scatter b_obs12 k, mcolor(navy) msymbol(circle) msize(vlarge)),          ///
    yline(0, lcolor(gs8) lpattern(dash))                                         ///
    xline(-0.5, lcolor(gs10) lpattern(shortdash))                                ///
    text(-0.22 -0.3 "Observed Y, both worlds", color(navy) size(large) placement(e)) ///
    text(-0.32 -0.3 "Y(0), World 2", color(maroon) size(large) placement(e))     ///
    text(-0.42 -0.3 "Y(0), World 1", color(gs6) size(large) placement(e))        ///
    xtitle("Event year (k = -1 omitted)", size(large))                           ///
    ytitle("Event-study coefficient", size(large))                               ///
    xlabel(-3(1)3, labsize(large)) xscale(range(-3.35 3.35))                                               ///
    ylabel(-0.4(0.2)0.8, angle(horizontal) format(%3.1f) labsize(large)) yscale(range(-0.48 0.8))         ///
    title("Worlds 1 and 2: the pre-trend test passes (p = `fpF_y_w1')", size(large))      ///
    subtitle("Seed 2024, 4,000 firms; 95 percent intervals", size(large))        ///
    legend(off) scheme(s1color) xsize(6) ysize(4.2)
graph export figures/fig1a_sv_parallel_trends_w12.png, replace width(2000)

* Figure 1b: World 3 (the test rejects; the PTA holds)
twoway (rcap lo_w3 hi_w3 k_r, lcolor(gs8))                                       ///
       (scatter b_w3 k_r, mcolor(gs6) msymbol(circle_hollow) msize(large))       ///
       (rcap lo_obs3 hi_obs3 k, lcolor(navy))                                    ///
       (scatter b_obs3 k, mcolor(navy) msymbol(circle) msize(vlarge)),           ///
    yline(0, lcolor(gs8) lpattern(dash))                                         ///
    xline(-0.5, lcolor(gs10) lpattern(shortdash))                                ///
    text(-0.24 -0.3 "Observed Y, World 3", color(navy) size(large) placement(e)) ///
    text(-0.35 -0.3 "Y(0), World 3", color(gs6) size(large) placement(e))        ///
    xtitle("Event year (k = -1 omitted)", size(large))                           ///
    ytitle("Event-study coefficient", size(large))                               ///
    xlabel(-3(1)3, labsize(large)) xscale(range(-3.35 3.35))                                               ///
    ylabel(-0.4(0.2)0.8, angle(horizontal) format(%3.1f) labsize(large)) yscale(range(-0.48 0.8))         ///
    title("World 3: the pre-trend test rejects (p < 0.001)", size(large))            ///
    subtitle("Seed 2024, 4,000 firms; 95 percent intervals", size(large))        ///
    legend(off) scheme(s1color) xsize(6) ysize(4.2)
graph export figures/fig1b_sv_parallel_trends_w3.png, replace width(2000)
restore

di "Done: tables/tab_sv_parallel_trends_worlds.tex, figures/fig1a_sv_parallel_trends_w12.png, figures/fig1b_sv_parallel_trends_w3.png"


********************************************************************************
* MODULE 2: the pre-trend test as a screen (Monte Carlo)
*
* What it demonstrates. Even when the violation of parallel trends is present
* before treatment, so that observed outcomes carry information about it, the
* joint pre-trend test may reject far less often than the DiD t-test rejects
* a true null of no effect. Conditioning on a passed pre-trend test does not
* remove the bias in the event-study coefficients and may enlarge it
* (Roth 2022).
*
* DGP. N = 200 firms (100 treated), event years k = -3,...,+3, treatment at
* k = 0. alpha_i ~ N(D_i, 1); lambda_k = 0.1*k; eps_ik ~ N(0, 1).
*   Y(0) = alpha + lambda + gamma*k*D + eps,  Y(1) = Y(0)  (ATT = 0).
* The untreated potential outcomes of treated firms drift by gamma per year
* relative to control firms in every year, before and after treatment.
* gamma in {0, 0.025, ..., 0.15}; 2,000 samples for each gamma.
*
* Analytical benchmarks (checked at tolerance 0.01):
*   E[static TWFE delta] = gamma*(1.5 - (-2)) = 3.5*gamma
*   E[beta_k]            = gamma*(k + 1), so the post average is 2.5*gamma
*
* Estimation: xtreg, fe with firm-clustered standard errors.
*
* Requirements: base Stata. Output: tables/tab_sv_parallel_trends_pretest.tex,
* figures/fig2_sv_parallel_trends_power.png
********************************************************************************
clear all
set seed 2024

tempname sims
tempfile mc
postfile `sims' gamma d_hat d_rej pre_p pre_rej avgpost avg_rej           ///
    b1 b2 b4 b5 b6 b7 using `mc', replace

foreach g in 0 0.025 0.05 0.075 0.10 0.125 0.15 {
    forvalues r = 1/2000 {
        quietly {
            clear
            set obs 200
            gen firm = _n
            gen D = firm <= 100
            gen alpha = D + rnormal()
            expand 7
            bysort firm: gen t = _n
            gen k = t - 4
            gen y = alpha + 0.1*k + `g'*k*D + rnormal()
            gen DPost = D*(k >= 0)
            xtset firm t
            forvalues j = 1/7 {
                gen dk`j' = D*(t==`j')
            }
            xtreg y i.t DPost, fe vce(cluster firm)
            local d  = _b[DPost]
            local dr = abs(_b[DPost]/_se[DPost]) > invttail(e(N_clust) - 1, 0.025)
            xtreg y i.t dk1 dk2 dk4 dk5 dk6 dk7, fe vce(cluster firm)
            testparm dk1 dk2
            local pp = r(p)
            lincom (dk4 + dk5 + dk6 + dk7)/4
            local ap  = r(estimate)
            local apr = r(p) < 0.05
            post `sims' (`g') (`d') (`dr') (`pp') (`pp' < 0.05) (`ap') (`apr') ///
                (_b[dk1]) (_b[dk2]) (_b[dk4]) (_b[dk5]) (_b[dk6]) (_b[dk7])
        }
    }
    di "gamma = `g' done"
}
postclose `sims'
use `mc', clear
gen pass = 1 - pre_rej

* ---- benchmarks ----
di _n "Monte Carlo benchmarks (tolerance 0.01):"
foreach g in 0 0.05 0.10 0.15 {
    quietly summarize d_hat if abs(gamma - `g') < 1e-6
    di "  gamma = " %5.3f `g' ": mean delta = " %6.3f r(mean) "  analytical 3.5*gamma = " %6.3f 3.5*`g'
    assert abs(r(mean) - 3.5*`g') < 0.01
    quietly summarize avgpost if abs(gamma - `g') < 1e-6
    di "               mean post avg = " %6.3f r(mean) "  analytical 2.5*gamma = " %6.3f 2.5*`g'
    assert abs(r(mean) - 2.5*`g') < 0.01
}

* ---- the table at gamma = 0.10 ----
di _n "Numbers for the slide (gamma = 0.10, true ATT = 0):"
quietly summarize pre_rej if abs(gamma - 0.10) < 1e-6
scalar s_prerej = r(mean)
quietly summarize d_rej if abs(gamma - 0.10) < 1e-6
scalar s_drej = r(mean)
quietly summarize d_rej if abs(gamma - 0.10) < 1e-6 & pass==1
scalar s_drej_p = r(mean)
quietly summarize d_hat if abs(gamma - 0.10) < 1e-6
scalar s_d = r(mean)
quietly summarize d_hat if abs(gamma - 0.10) < 1e-6 & pass==1
scalar s_d_p = r(mean)
quietly summarize avgpost if abs(gamma - 0.10) < 1e-6
scalar s_ap = r(mean)
quietly summarize avgpost if abs(gamma - 0.10) < 1e-6 & pass==1
scalar s_ap_p = r(mean)
quietly summarize avg_rej if abs(gamma - 0.10) < 1e-6
scalar s_aprej = r(mean)
quietly summarize avg_rej if abs(gamma - 0.10) < 1e-6 & pass==1
scalar s_aprej_p = r(mean)
di "  pre-trend test rejects                       : " %5.3f s_prerej
di "  static DiD t-test rejects ATT = 0, all       : " %5.3f s_drej
di "  static DiD t-test rejects ATT = 0, passers   : " %5.3f s_drej_p
di "  mean static delta, all / passers (true 0.35) : " %6.4f s_d " / " %6.4f s_d_p
di "  mean post avg, all / passers (true 0.25)     : " %5.3f s_ap " / " %5.3f s_ap_p
di "  post avg test rejects, all / passers         : " %5.3f s_aprej " / " %5.3f s_aprej_p
quietly summarize pre_rej if gamma == 0
di "  size check at gamma = 0: pre-trend test      : " %5.3f r(mean)
quietly summarize d_rej if gamma == 0
di "  size check at gamma = 0: static DiD t-test   : " %5.3f r(mean)

foreach s in s_prerej s_drej s_drej_p s_d s_d_p s_ap s_ap_p s_aprej s_aprej_p {
    local f`s' : di %4.2f `s'
    local f`s' = strtrim("`f`s''")
}
file open tab using tables/tab_sv_parallel_trends_pretest.tex, write replace
file write tab "\begin{tabular}{@{}lcc@{}}" _n
file write tab "\toprule" _n
file write tab " & All samples & Pre-test passed \\" _n
file write tab "\midrule" _n
file write tab "Share of samples & 1.00 & \textcolor{Maroon}{\bfseries " %4.2f (1 - s_prerej) "} \\" _n
file write tab "Mean \$\hat\delta\$ (ATT \$= 0\$) & `fs_d' & `fs_d_p' \\" _n
file write tab "\$t\$-test rejects ATT \$= 0\$ & `fs_drej' & \textcolor{Maroon}{\bfseries `fs_drej_p'} \\" _n
file write tab "Mean post \$\hat\beta_k\$ (bias 0.25) & `fs_ap' & \textcolor{Maroon}{\bfseries `fs_ap_p'} \\" _n
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_sv_parallel_trends_pretest.tex

* ---- figure 2: rejection rates against the size of the violation ----
preserve
collapse (mean) pre_rej d_rej avg_rej, by(gamma)
list, clean noobs
twoway (connected d_rej gamma, lcolor(maroon) mcolor(maroon) lwidth(thick) msymbol(diamond))   ///
       (connected pre_rej gamma, lcolor(navy) mcolor(navy) lwidth(thick) msymbol(circle)),     ///
    yline(0.05, lcolor(gs8) lpattern(dash))                                                    ///
    text(0.93 0.072 "DiD t-test rejects a true" "ATT of zero", color(maroon) size(large) placement(w) justification(right)) ///
    text(0.33 0.149 "Pre-trend F test rejects" "parallel pre-trends", color(navy) size(large) placement(w) justification(right)) ///
    xtitle("Drift in Y(0) per year (gamma), in s.d. of the error", size(medlarge))                ///
    ytitle("Share of 2,000 samples", size(medlarge))                                                           ///
    xlabel(0(0.05)0.15, format(%3.2f) labsize(medlarge)) xscale(range(0 0.157)) ylabel(0(0.2)1, angle(horizontal) format(%3.1f) labsize(medlarge))       ///
    title("The pre-trend test is the weaker of the two tests", size(large))                                 ///
    subtitle("Seed 2024: 200 firms, true ATT = 0; dashed line = 0.05", size(medlarge))            ///
    legend(off) scheme(s1color) xsize(7) ysize(4.6)
graph export figures/fig2_sv_parallel_trends_power.png, replace width(2000)
restore

di "Done: tables/tab_sv_parallel_trends_pretest.tex, figures/fig2_sv_parallel_trends_power.png"

log close
