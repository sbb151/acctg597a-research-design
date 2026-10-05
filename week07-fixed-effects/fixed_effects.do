********************************************************************************
* ACCTG 597A -- Week 7: Fixed Effects, the Variation They Leave and the
* Comparison They Draw
* Mummolo, J., and E. Peterson (2018), "Improving the Interpretation of Fixed
*   Effects Regression Results," Political Science Research and Methods 6(4): 829-835.
* Kropko, J., and R. Kubinec (2020), "Interpretation and Identification of
*   Within-Unit and Cross-Sectional Variation in Panel Data Models," PLOS ONE
*   15(4): e0231349.
* Running example: Anantharaman, D., and Y. G. Lee (2014), "Managerial Risk
*   Taking Incentives and Corporate Pension Policy," Journal of Financial
*   Economics 111(2): 328-351.
*
* Self-contained simulation, base Stata only (version 17 or later). Run from
* this folder:  do fixed_effects.do
* Figures are written to figures/ and a booktabs table fragment to tables/.
* Runtime under a minute.
*
* Every module simulates the same firm-year panel of CFO vega (900 firms x 5
* years, seed 2014; 25 percent of vega's sum of squares within firms, an
* assumption, because A&L report no within statistics; 72 firms whose vega
* never changes), standardized by its pooled mean and SD as A&L do.
*   1. Within variation (Mummolo and Peterson). Underfunding responds to vega
*      with slope 0.0176 (A&L Table 4A, Model 3, at the 95th percentile of
*      distress). Dummies, xtreg, and the regression on residualized vega give
*      one slope; SST = SSW + SSB; the one-way slope is a weighted average of
*      firm slopes with weight zero for constant-vega firms; SD(x~) is 0.50
*      (firm effects) and 0.47 (firm and year effects) of the pooled SD.
*   2. Which comparison (Kropko and Kubinec). Vega splits exactly into a
*      between-firm, a year, and a residual part, each with its own slope in
*      the pension equity share. Pooled, year, firm, and two-way slopes equal
*      share-weighted averages of the kept parts' slopes; firm and year FE
*      equal weighted averages of subset slopes (K&K Eqs 8, 10-11); two-way FE
*      equals K&K Eqs 12 and 13; firm FE minus two-way FE is a year term.
*   3. Identification. The sampling SD of the two-way slope as the residual
*      part of vega shrinks (Monte Carlo, 500 draws per point).
*   4. K&K's homogeneous design: x is additive, Stata drops a dummy, and the
*      reported slope depends on the order of the dummies.
* Identity checks print and assert; each should be zero to machine precision.
********************************************************************************

version 17
clear all
set more off
capture log close
log using fixed_effects.log, replace text
capture mkdir figures
capture mkdir tables

* The shared panel, built after each set seed
capture program drop fe_panel
program define fe_panel
    local N      900
    local T      5
    local pconst 0.08
    local target 0.25

    clear
    set obs `N'
    gen int firm = _n
    gen double mu = rnormal(0, 1)
    gen byte constant = (_n <= round(`pconst' * `N'))
    expand `T'
    bysort firm: gen int year = _n
    gen double lam = 0.25*(year==1) + 0.20*(year==2) + 0.10*(year==3) ///
                   - 0.25*(year==4) - 0.30*(year==5)
    gen double nu = rnormal(0, 1)
    replace nu = 0 if constant
    gen double A = mu + lam
    replace A = -1.6 if constant

    * Solve SSW(k) = target * SST(k) for the scale k on nu (a quadratic in k)
    foreach v in A nu {
        bysort firm: egen double `v'_fm = mean(`v')
        quietly summarize `v', meanonly
        gen double `v'_w = `v' - `v'_fm
        gen double `v'_t = `v' - r(mean)
    }
    foreach p in w t {
        gen double aa_`p' = A_`p'^2
        gen double an_`p' = A_`p' * nu_`p'
        gen double nn_`p' = nu_`p'^2
        foreach s in aa an nn {
            quietly summarize `s'_`p', meanonly
            scalar Z`s'`p' = r(sum)
        }
    }
    * (scalars carry a Z prefix so that no variable name shadows them)
    scalar qa = Znnw - `target' * Znnt
    scalar qb = 2 * (Zanw - `target' * Zant)
    scalar qc = Zaaw - `target' * Zaat
    scalar k_scale = (-qb + sqrt(qb^2 - 4*qa*qc)) / (2*qa)
    gen double vraw = A + k_scale * nu
    drop A_* nu_* aa_* an_* nn_* A

    quietly summarize vraw
    gen double vega = (vraw - r(mean)) / r(sd)
    label variable vega "CFO vega, standardized by pooled mean and SD"

    * Realized within share (must equal the target)
    bysort firm: egen double v_fm = mean(vega)
    gen double v_w2 = (vega - v_fm)^2
    gen double v_t2 = vega^2
    quietly summarize v_w2, meanonly
    scalar ssw_v = r(sum)
    quietly summarize v_t2, meanonly
    scalar sst_v = r(sum)
    scalar share_w = ssw_v / sst_v
    drop v_fm v_w2 v_t2
    di "Panel: `N' firms x `T' years; scale k on nu = " %8.6f k_scale
    di "Within-firm share of the sum of squares of vega (target `target'): " %10.8f share_w
    assert abs(share_w - `target') < 1e-8
    xtset firm year
end


********************************************************************************
* MODULE 1: the variation a fixed effects coefficient uses (Mummolo and Peterson)
********************************************************************************
clear
set seed 2014
fe_panel

* ---- outcome ----------------------------------------------------------------
bysort firm (year): gen double a_i = -0.012*mu + rnormal(0, 0.15) if _n == 1
bysort firm (year): replace a_i = a_i[1]
gen double uf = 0.158 + a_i + 0.0176*vega + rnormal(0, 0.20)
label variable uf "Underfunding (fraction of obligations)"

* ---- (1) three routes to the firm-and-year slope ------------------------------
quietly regress uf vega i.firm i.year
scalar b_dum = _b[vega]
quietly xtreg uf vega i.year, fe vce(cluster firm)
scalar b_xt  = _b[vega]
scalar t_xt  = _b[vega] / _se[vega]
quietly regress vega i.firm i.year
predict double vt2, residuals              // x~ under firm and year effects
quietly regress uf vt2
scalar b_fwl = _b[vt2]
di _n "Identity checks (should be zero to machine precision):"
di "  (1a) dummies minus xtreg fe        : " %12.10f (b_dum - b_xt)
di "  (1b) dummies minus FWL residual    : " %12.10f (b_dum - b_fwl)
assert abs(b_dum - b_xt) < 1e-6
assert abs(b_dum - b_fwl) < 1e-6

* ---- (2) one-way within slope as a ratio of within sums (M&P Eq. 2) ----------
bysort firm: egen double v_fm = mean(vega)
bysort firm: egen double u_fm = mean(uf)
gen double vw = vega - v_fm                // x~ under firm effects only
gen double uw = uf - u_fm
gen double xy = vw*uw
gen double xx = vw^2
quietly summarize xy, meanonly
scalar sxy = r(sum)
quietly summarize xx, meanonly
scalar sxx = r(sum)
quietly xtreg uf vega, fe
scalar b_fe1 = _b[vega]
di "  (2)  xtreg fe minus sum(xw*yw)/sum(xw^2): " %12.10f (b_fe1 - sxy/sxx)
assert abs(b_fe1 - sxy/sxx) < 1e-6

* ---- (3) SST = SSW + SSB -----------------------------------------------------
quietly summarize vega
scalar sst = r(Var) * (r(N) - 1)
scalar nobs = r(N)
gen double bb = (v_fm - r(mean))^2
quietly summarize bb, meanonly
scalar ssb = r(sum)
di "  (3)  SST - SSW - SSB               : " %12.10f (sst - sxx - ssb)
assert abs(sst - sxx - ssb) < 1e-6

* ---- (4) firm weights --------------------------------------------------------
preserve
collapse (sum) sxy_i = xy sxx_i = xx (max) constant, by(firm)
gen double b_i = sxy_i / sxx_i               // missing for constant firms
quietly summarize sxx_i if constant, meanonly
scalar w_const = r(max)
gen double wb = sxx_i * b_i
quietly summarize wb, meanonly
scalar num_w = r(sum)
quietly summarize sxx_i, meanonly
scalar den_w = r(sum)
quietly count if constant
scalar n_const = r(N)
quietly count
scalar n_firms = r(N)
restore
di "  (4a) b_FE minus sum(w_i b_i)/sum(w_i): " %12.10f (b_fe1 - num_w/den_w)
di "  (4b) largest weight of a constant-vega firm: " %12.10f w_const
assert abs(b_fe1 - num_w/den_w) < 1e-6
assert w_const == 0

* ---- (5) standard deviations and described effects ---------------------------
quietly summarize vega
scalar sd_x  = r(sd)
quietly summarize vw
scalar sd_x1 = r(sd)
quietly summarize vt2
scalar sd_x2 = r(sd)
quietly summarize vt2
scalar ss2 = r(Var) * (r(N) - 1)
scalar share_2 = ss2 / sst
di "  (5a) SD(x~1)/SD(x) minus sqrt(SSW/SST): " %12.10f (sd_x1/sd_x - sqrt(sxx/sst))
di "  (5b) SD(x~2)/SD(x) minus sqrt(SS2/SST): " %12.10f (sd_x2/sd_x - sqrt(share_2))
* With year effects as well, a constant-vega firm's residual is -(xbar_t - xbar), not zero:
* it stays in the estimation as a comparison firm. Its share of the two-way sum of squares:
gen double vt2sq = vt2^2
quietly summarize vt2sq if constant, meanonly
scalar ss2_const = r(sum)
quietly summarize vt2sq, meanonly
scalar sh2_const = ss2_const / r(sum)
di "  constant-vega firms' share of the two-way sum of squares: " %8.6f sh2_const
assert sh2_const > 0      // with year effects, constant firms remain as comparison firms
assert abs(sd_x1/sd_x - sqrt(sxx/sst)) < 1e-6
assert abs(sd_x2/sd_x - sqrt(share_2)) < 1e-6

* cross-sectional comparison (year effects only), for contrast
quietly regress uf vega i.year, vce(cluster firm)
scalar b_cs = _b[vega]
scalar t_cs = _b[vega] / _se[vega]

* A&L's 50th-to-95th percentile shift in CFO vega: 1.993 pooled SDs (Table 1)
scalar al_shift = (202.77 - 19.87) / 91.77
scalar shift_in_sd2 = al_shift / (sd_x2 / sd_x)

di _n "Numbers for the slides:"
quietly summarize uf
di "  SD of underfunding (A&L Table 1: 0.247): " %6.4f r(sd)
di "  N firm-years                      : " nobs
di "  firms with constant vega          : " n_const " of " n_firms
di "  within-firm share of SST (one-way): " %6.4f (sxx/sst)
di "  residual share, firm and year FE  : " %6.4f share_2
di "  SD(x) pooled                      : " %6.4f sd_x
di "  SD(x~) firm effects               : " %6.4f sd_x1
di "  SD(x~) firm and year effects      : " %6.4f sd_x2
di "  slope, firm and year FE           : " %7.5f b_xt "  (t = " %4.2f t_xt ")"
di "  slope, year FE only               : " %7.5f b_cs "  (t = " %4.2f t_cs ")"
di "  slope, firm FE only (one-way)     : " %7.5f b_fe1
di "  A&L Table 4A M3 slope at distress P95 (x 1,000): " %6.2f (22.0 + 2.340*(-1.877))
di "  described effect, pooled SD       : " %7.5f (b_xt*sd_x)
di "  described effect, SD(x~) two-way  : " %7.5f (b_xt*sd_x2)
di "  A&L 50th-95th shift, pooled SDs   : " %6.4f al_shift
di "  same shift in two-way SD(x~) units: " %6.4f shift_in_sd2


* Figure: the distribution of vega before and after firm and year effects
* (the analogue of Mummolo and Peterson's Figure 1, left panel)
local s1 : di %4.2f sd_x
local s2 : di %4.2f sd_x2
twoway (kdensity vega, lcolor(maroon) lwidth(medthick))                     ///
       (kdensity vt2, lcolor(navy) lwidth(thick)),                         ///
    text(0.30 2.2 "pooled, SD = `s1'", color(maroon) size(medium) placement(e)) ///
    text(0.80 0.9 "after firm and year effects, SD = `s2'", color(navy) size(medium) placement(e)) ///
    xtitle("CFO vega, pooled standard deviations") ytitle("Density")       ///
    xlabel(-3(1)4) ylabel(, angle(horizontal))                             ///
    title("Firm and year effects leave half the spread")                   ///
    subtitle("Simulated panel (seed 2014), 900 firms x 5 years")           ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig1_within_variation.png, replace width(2000)


********************************************************************************
* MODULE 2: which comparison each fixed effects choice draws (Kropko and Kubinec)
********************************************************************************
clear
set seed 2014
fe_panel

local theta 0.0183
local phi   0.03
local gamma 0.0026

* ---- components of vega --------------------------------------------------------
quietly summarize vega, meanonly
scalar xbar = r(mean)
bysort firm: egen double xb_i = mean(vega)
bysort year: egen double xb_t = mean(vega)
gen double B  = xb_i - xbar
gen double Yc = xb_t - xbar
gen double R  = vega - xb_i - xb_t + xbar
foreach c in B Yc R {
    gen double `c'2 = `c'^2
    quietly summarize `c'2, meanonly
    scalar SS`c' = r(sum)
}
scalar SStot = SSB + SSYc + SSR
di "Shares of vega's sum of squares: between " %6.4f (SSB/SStot) ///
   ", year " %6.4f (SSYc/SStot) ", two-way residual " %6.4f (SSR/SStot)

* ---- outcome -------------------------------------------------------------------
bysort firm (year): gen double a_i = rnormal(0, 0.10) if _n == 1
bysort firm (year): replace a_i = a_i[1]
gen double sys = `theta'*B + `phi'*Yc + `gamma'*R        // systematic part
gen double eq  = 0.603 + a_i + sys + rnormal(0, 0.10)
label variable eq "Pension equity share"

* ---- estimators -----------------------------------------------------------------
program define fe4, rclass
    * returns pooled, year FE, firm FE, two-way FE slopes of `1' on vega
    quietly regress `1' vega
    return scalar pool = _b[vega]
    quietly regress `1' vega i.year
    return scalar yfe = _b[vega]
    quietly xtreg `1' vega, fe
    return scalar ffe = _b[vega]
    quietly xtreg `1' vega i.year, fe
    return scalar tw = _b[vega]
end
fe4 sys
foreach e in pool yfe ffe tw {
    scalar S_`e' = r(`e')
}
fe4 eq
foreach e in pool yfe ffe tw {
    scalar E_`e' = r(`e')
}
quietly regress eq vega i.year, vce(cluster firm)
scalar t_yfe = _b[vega] / _se[vega]
quietly xtreg eq vega i.year, fe vce(cluster firm)
scalar t_tw = _b[vega] / _se[vega]

di _n "Identity checks (should be zero to machine precision):"
* (4) systematic part: each estimator = SS-weighted average of kept slopes
scalar L_pool = (`theta'*SSB + `phi'*SSYc + `gamma'*SSR) / SStot
scalar L_yfe  = (`theta'*SSB + `gamma'*SSR) / (SSB + SSR)
scalar L_ffe  = (`phi'*SSYc + `gamma'*SSR) / (SSYc + SSR)
scalar L_tw   = `gamma'
foreach e in pool yfe ffe tw {
    di "  (4) `e' on systematic part minus weighted average of kept slopes: " %12.10f (S_`e' - L_`e')
    assert abs(S_`e' - L_`e') < 1e-6
}

* ---- (1) firm FE as a weighted average of firm slopes (K&K Eqs 9-11) --------------
bysort firm: egen double y_i = mean(eq)
gen double xw = vega - xb_i
gen double yw = eq - y_i
gen double xwyw = xw*yw
gen double xw2  = xw^2
* ---- (2) year FE as a weighted average of year slopes (K&K Eq 8) ------------------
bysort year: egen double y_t = mean(eq)
gen double xy_ = vega - xb_t
gen double yy_ = eq - y_t
gen double xyyy = xy_*yy_
gen double xy2  = xy_^2
* ---- (3) two-way: Eq 12 and Eq 13 ---------------------------------------------------
quietly summarize eq, meanonly
scalar ybar = r(mean)
gen double ytw = eq - y_i - y_t + ybar
gen double Rytw = R*ytw
preserve
    collapse (sum) n1 = xwyw d1 = xw2, by(firm)
    gen double b_i = n1/d1
    gen double wb = d1*b_i
    quietly summarize wb, meanonly
    scalar num1 = r(sum)
    quietly summarize d1, meanonly
    scalar den1 = r(sum)
restore
preserve
    collapse (sum) n2 = xyyy d2 = xy2, by(year)
    gen double b_t = n2/d2
    gen double wb = d2*b_t
    quietly summarize wb, meanonly
    scalar num2 = r(sum)
    quietly summarize d2, meanonly
    scalar den2 = r(sum)
restore
quietly summarize Rytw, meanonly
scalar n12 = r(sum)
* Eq 13: within each year, regress firm-demeaned y on firm-demeaned x
gen double ys = eq - y_i
bysort year: egen double ys_t = mean(ys)
bysort year: egen double xs_t = mean(xw)
gen double xs_c = xw - xs_t
gen double ys_c = ys - ys_t
gen double n13 = xs_c*ys_c
gen double d13 = xs_c^2
preserve
    collapse (sum) n13 d13, by(year)
    gen double b13 = n13/d13
    gen double wb = d13*b13
    quietly summarize wb, meanonly
    scalar num13 = r(sum)
    quietly summarize d13, meanonly
    scalar den13 = r(sum)
restore
di "  (1) firm FE minus sum(w_i b_i)/sum(w_i)      : " %12.10f (E_ffe - num1/den1)
di "  (2) year FE minus sum(w_t b_t)/sum(w_t)      : " %12.10f (E_yfe - num2/den2)
di "  (3a) two-way FE minus K&K Eq 12 ratio         : " %12.10f (E_tw - n12/SSR)
di "  (3b) two-way FE minus K&K Eq 13 average       : " %12.10f (E_tw - num13/den13)
assert abs(E_ffe - num1/den1) < 1e-6
assert abs(E_yfe - num2/den2) < 1e-6
assert abs(E_tw - n12/SSR) < 1e-6
assert abs(E_tw - num13/den13) < 1e-6

* ---- (5) the gap between firm FE and two-way FE (appendix B1) ------------------------
quietly regress eq vega i.firm i.year
gen double d_hat = 0
forvalues t = 2/5 {
    quietly replace d_hat = _b[`t'.year] if year == `t'
}
gen double gap_term = d_hat*Yc
quietly summarize gap_term, meanonly
scalar gap_rhs = r(sum) / (SSYc + SSR)
di "  (5) (firm FE - two-way FE) minus N sum d_t (xbar_t - xbar)/SSW: " %12.10f ((E_ffe - E_tw) - gap_rhs)
assert abs((E_ffe - E_tw) - gap_rhs) < 1e-6

di _n "Numbers for the slides (x 1,000, as A&L report them):"
di "  shares: between " %6.4f (SSB/SStot) "  year " %6.4f (SSYc/SStot) "  residual " %6.4f (SSR/SStot)
di "  systematic part: pooled " %6.2f (1000*S_pool) "  year FE " %6.2f (1000*S_yfe) ///
   "  firm FE " %6.2f (1000*S_ffe) "  two-way " %6.2f (1000*S_tw)
di "  simulated eq   : pooled " %6.2f (1000*E_pool) "  year FE " %6.2f (1000*E_yfe) ///
   "  firm FE " %6.2f (1000*E_ffe) "  two-way " %6.2f (1000*E_tw)
di "  t (clustered by firm): year FE " %5.2f t_yfe "  two-way " %5.2f t_tw
quietly summarize eq
di "  SD of equity share (A&L Table 1: 0.152): " %6.4f r(sd)
di "  A&L Model 3 slopes at distress P95 (x 1,000): Table 5 " %6.2f (18.6 + 2.050*(-1.877)) "  Table 6A " %6.2f (4.510 + 1.010*(-1.877))


* Table: the four fixed effects choices
file open tab using tables/tab_fe_slices.tex, write replace
file write tab "\begin{tabular}{@{}lccc@{}}" _n
file write tab "\toprule" _n
file write tab "Fixed effects & Kept parts of vega & Slope, systematic part & Estimate \\" _n
file write tab "\midrule" _n
foreach e in pool yfe ffe tw {
    local lab = cond("`e'"=="pool","None",cond("`e'"=="yfe","Year",cond("`e'"=="ffe","Firm","Firm and year")))
    local kept = cond("`e'"=="pool","all",cond("`e'"=="yfe","between, residual",cond("`e'"=="ffe","year, residual","residual")))
    local a : di %4.1f (1000*S_`e')
    local b : di %4.1f (1000*E_`e')
    file write tab "`lab' & `kept' & `=strtrim("`a'")' & `=strtrim("`b'")' \\" _n
}
file write tab "\bottomrule" _n
file write tab "\end{tabular}" _n
file close tab
type tables/tab_fe_slices.tex


********************************************************************************
* MODULE 3: the two-way slope as the residual part of vega shrinks
********************************************************************************
* ---- Monte Carlo: the two-way slope as the residual component shrinks ------------
mata:
void mc_twfe(real scalar N, real scalar T, real scalar R, real rowvector cs,
             real scalar gam, real scalar se)
{
    real scalar j, r, c, ss, b
    real matrix X, Y, Xt, Yt, out
    real colvector mu, a
    real rowvector lam, bs, an, sh
    lam = (0.25, 0.20, 0.10, -0.25, -0.30)
    out = J(cols(cs), 5, .)
    for (j = 1; j <= cols(cs); j++) {
        c = cs[j]
        bs = J(1, R, .); an = J(1, R, .); sh = J(1, R, .)
        for (r = 1; r <= R; r++) {
            mu = rnormal(N, 1, 0, 1)
            a  = rnormal(N, 1, 0, 0.10)
            X  = (mu :+ c :* rnormal(N, T, 0, 1)) :+ lam   // column first, then row (Mata colon rules)
            X  = (X :- mean(vec(X))) :/ sqrt(variance(vec(X)))
            Y  = a :+ gam :* X :+ rnormal(N, T, 0, se)
            Xt = X :- rowsum(X)/T :- colsum(X)/N :+ sum(X)/(N*T)
            Yt = Y :- rowsum(Y)/T :- colsum(Y)/N :+ sum(Y)/(N*T)
            ss = sum(Xt:^2)
            bs[r] = sum(Xt:*Yt)/ss
            an[r] = se/sqrt(ss)
            sh[r] = ss/((N*T) - 1)
        }
        out[j, .] = (c, mean(sh'), mean(bs'), sqrt(variance(bs')), mean(an'))
    }
    st_matrix("MC", out)
}
end
mata: mc_twfe(900, 5, 500, (0.70, 0.35, 0.18, 0.09, 0.045, 0.022, 0.011, 0.0055), `gamma', 0.10)
matrix colnames MC = scale share mean_b sd_b analytic
matrix list MC, format(%9.5f)
di _n "Monte Carlo check: simulated SD of the two-way slope versus sigma_e/sqrt(sum R^2)"
forvalues j = 1/`=rowsof(MC)' {
    scalar ratio = MC[`j', 4] / MC[`j', 5]
    di "  residual share " %7.5f MC[`j', 2] ": SD " %8.5f MC[`j', 4] "  analytic " %8.5f MC[`j', 5] "  ratio " %5.3f ratio
    assert abs(ratio - 1) < 0.10
}


* Figure: sampling SD of the two-way slope against the residual share
clear
svmat MC, names(col)
gen sd1000 = 1000*sd_b
gen an1000 = 1000*analytic
twoway (line an1000 share, lcolor(gs8))                                    ///
       (scatter sd1000 share, mcolor(navy) msymbol(circle) msize(large)), ///
    xscale(log reverse) yscale(log)                                        ///
    xlabel(0.1 "10%" 0.01 "1%" 0.001 "0.1%" 0.0001 "0.01%")                ///
    ylabel(3 10 30 100 300, angle(horizontal))                             ///
    yline(`=1000*`gamma'', lcolor(maroon) lpattern(dash))                  ///
    xtitle("Share of vega's variance that is neither firm nor year")       ///
    ytitle("SD of the two-way slope (x 1,000)")                            ///
    title("Less residual variation, a noisier two-way slope")              ///
    subtitle("500 draws per point; dashed line: the true slope")           ///
    legend(off) scheme(s1color) xsize(9) ysize(5)
graph export figures/fig2_twoway_precision.png, replace width(2000)


********************************************************************************
* MODULE 4: K&K's homogeneous design, where the two-way slope is not identified
********************************************************************************
* ---- K&K's DGP with homogeneous slopes (Eqs 20-21): not identified -------------------
preserve
    clear
    set obs 30
    gen int case = _n
    gen double alpha_i = rnormal()
    expand 30
    bysort case: gen int time = _n
    gen double alpha_t = rnormal() if case == 1
    bysort time (case): replace alpha_t = alpha_t[1]
    local bta 1
    local gma -3
    gen double x = (alpha_i - alpha_t) / (`bta' - `gma')
    gen double y = (`bta'*alpha_i - `gma'*alpha_t) / (`bta' - `gma')   // noise-free, as in K&K Eq 21
    regress y x i.case i.time
    scalar kk_b1 = _b[x]
    regress y x i.time i.case
    scalar kk_b2 = _b[x]
    di "  K&K homogeneous DGP (beta_t = `bta', gamma_i = `gma'):"
    di "    order x i.case i.time -> b = " %8.4f kk_b1
    di "    order x i.time i.case -> b = " %8.4f kk_b2
restore


di _n "Done: fixed_effects.do"
log close
