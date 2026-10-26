*==============================================================================
* shift_share_simulation.do: a shift-share ("simulated instrument") measure of the
* scope of size-based reporting regulation, as in Breuer, Leuz, and
* Vanhaverbeke (2025, JAE), "Reporting Regulation and Corporate Innovation".
* Teaching note, ACCTG 597A Week 10.
*
* The design in one line. An industry-country-year cell (k, c, t) has a share
* x of firms above the country-year size threshold tau_ct that report (the
* scope of the mandate). The outcome y is innovation in the cell. Firm size,
* and so x, responds to the same innovation shocks that move y (reverse
* causality), so OLS of y on x is biased. The shift-share measure freezes the
* size distribution: it is the share of the industry's BASE-PERIOD firms that
* would lie above tau_ct,
*     z_kct = sum_b s_kb * g_bct,
* with s_kb the industry's pre-determined share of firms in size bin b (the
* "shares") and g_bct = 1[bin b lies above tau_ct] (the "shifts"). Only the
* thresholds move z; the firms' responses do not. Breuer, Leuz, and
* Vanhaverbeke regress y on z directly (a reduced form) with country-year and
* industry-year fixed effects, so that the identifying variation is the
* interaction of a country's threshold level with an industry's size
* distribution; the script also reports OLS on the actual share x and the 2SLS
* estimate that treats z as an instrument for x, which the paper does not run.
*
* DGP (seed 2025).
*   K = 30 industries, C = 10 countries, T = 8 years, n = 50 firms per cell;
*   base log size l_i0 ~ N(mu_k, 1), mu_k evenly spaced on [-1.2, 1.2];
*   six size bins cut at -1, -0.5, 0, 0.5, 1; thresholds tau_ct take these
*   values; countries cycle through the five levels in random order, so that
*   levels are distinct, and `nchange' of them move by +-0.5 in a random year
*   between 3 and 6 (three by default; the paper calls its thresholds
*   persistent, fn. 15);
*   current log size l_ict = l_i0 + lambda u_kct + nu_ict, nu ~ N(0, 0.3^2);
*   a firm above the threshold complies with the mandate with probability e;
*   x_kct = share of the cell's firms with l_ict >= tau_ct that comply;
*   y_kct = alpha_k + delta_ct + beta x_kct + u_kct
*           + gamma (w_kc - mean w) + rho_c small_k + e_kct,
*   u ~ N(0, 0.5^2) the innovation shock, e ~ N(0, 0.2^2); w_kc the share of
*   industry k's base-period firms in the size bin just above country c's
*   initial threshold, i.e., the incidence of a bundled size-based rule that
*   shares the thresholds (an audit exemption, a simplified format) and is not
*   the reporting mandate, a time-invariant country-industry factor correlated
*   with z but not identical to it; small_k the base share of firms below log
*   size 0; rho_c ~ N(0, sig_rho^2) a country-level factor whose incidence on
*   industries depends on their size structure.
*   lambda: reverse causality (innovation shocks grow firms into the mandate);
*   gamma:  an omitted country-industry factor correlated with the intensity
*           measure, the case the paper's identifying assumption rules out
*           (p. 9), other than through x;
*   sig_rho: exposure-correlated country factors, which average out only over
*           many countries (Borusyak, Hull, and Jaravel 2022) and are missed by
*           standard errors clustered below the country level;
*   e:      compliance with the mandate, which scales the first stage of x on z
*           (the paper's fn. 9: the estimate reflects average enforcement).
* Estimators: OLS of y on x, the reduced form of y on z, and 2SLS of y on x
*   instrumented by z, under two fixed-effect structures: (1) industry-year and
*   country-year effects, as in the paper, and (2) the same plus
*   country-industry effects, which leave only the threshold changes. Both are
*   exact demeanings on the balanced grid (two-way within each year; the
*   three-factor-interaction residual). Standard errors clustered by
*   country-year, two-way by country-industry and country-year (the paper's
*   scheme), and by country (the level of the shock); the first-stage F;
*   Rotemberg weights by size bin (Goldsmith-Pinkham, Sorkin, and Swift 2020)
*   and the analogous decomposition by country.
* Identities (asserted, one draw): ivregress and regress agree with the
*   demeaned ratios and with the Mata clustered standard error; 2SLS equals the
*   reduced form over the first stage; the weights sum to one and the weighted
*   bin (country) estimates sum to the 2SLS estimate; under structure (2) a z
*   built from the initial thresholds alone is demeaned to zero.
* Monte Carlo grids (R draws each, structure 1 with structure 2 alongside):
*   lambda, gamma, the number of countries, and e.
* Output: output/sim-bartik.js (the deck's data file, window.SIMBARTIK), four
*   figures in figures/ (twoway), and the log.
* Requirements: base Stata (Mata for the loops). About six minutes.
*==============================================================================
version 17
clear all
set more off
set seed 2025
local beta   -0.50
local R      300
local outjs  "output/sim-bartik.js"
capture mkdir figures
capture mkdir output

mata:
mata clear
mata set matastrict off

// ---------- exact demeaning on the balanced grid --------------------------------
// Cells are ordered industry-major: i = (k-1) C T + (c-1) T + t, so rowshape(v', K)
// is the K x (C T) matrix with column (c-1) T + t.
//   mode 1: industry-year and country-year effects (the paper): within each year,
//           subtract industry means (rows) and country means (columns) of the K x C grid.
//   mode 2: industry-year, country-year, and country-industry effects: the
//           three-factor-interaction residual
//           v - m_kc - m_kt - m_ct + m_k + m_c + m_t - m.
real colvector dm(real colvector v, real scalar K, real scalar C, real scalar T, real scalar mode)
{
    M  = rowshape(v', K)
    CT = C*T
    if (mode == 1) {
        for (t=1; t<=T; t++) {
            cols = range(0, C-1, 1)' :* T :+ t
            G = M[., cols]
            M[., cols] = G :- mean(G')' :- mean(G) :+ mean(vec(G))
        }
        return(vec(M'))
    }
    Mkc = J(K, CT, .); Mkt = J(K, CT, .); mc_ = J(K, CT, .); mt_ = J(K, CT, .)
    for (c=1; c<=C; c++) {
        cols = (c-1)*T :+ range(1, T, 1)'
        G = M[., cols]
        Mkc[., cols] = mean(G')' * J(1, T, 1)
        mc_[., cols] = J(K, T, mean(vec(G)))
    }
    for (t=1; t<=T; t++) {
        cols = range(0, C-1, 1)' :* T :+ t
        G = M[., cols]
        Mkt[., cols] = mean(G')' * J(1, C, 1)
        mt_[., cols] = J(K, C, mean(vec(G)))
    }
    Mct = J(K, 1, 1) * mean(M)
    mk  = mean(M')' * J(1, CT, 1)
    Rm  = M - Mkc - Mkt - Mct + mk + mc_ + mt_ :- mean(vec(M))
    return(vec(Rm'))
}

// ---------- one draw of the economy --------------------------------------------
void draw(real scalar K, real scalar C, real scalar T, real scalar n,
          real scalar lambda, real scalar gamma, real scalar sig_rho,
          real scalar nchange, real scalar e, real scalar beta,
          real colvector x, real colvector z, real colvector zlev, real colvector y,
          real matrix Zb, real colvector ctry, real colvector ind, real colvector cy,
          real colvector kc, real colvector iy, real matrix S, real colvector small,
          real matrix tau, real matrix Wb)
{
    grid  = (-1, -.5, 0, .5, 1)
    edges = (-1e9, grid)
    B     = cols(edges)
    mu    = rangen(-1.2, 1.2, K)
    sig_u = 0.5; sig_e = 0.2; sig_nu = 0.3
    // base-period sizes: one row per industry-country cell, n firms per row
    N0 = K*C
    L0 = J(N0, n, .)
    kk = J(N0, 1, .)
    for (k=1; k<=K; k++) {
        for (c=1; c<=C; c++) {
            i0 = (k-1)*C + c
            L0[i0, .] = mu[k] :+ rnormal(1, n, 0, 1)
            kk[i0] = k
        }
    }
    // base shares by industry, pooled across countries (pre-determined)
    S = J(K, B, 0)
    for (k=1; k<=K; k++) {
        Lk = vec(L0[select(1::N0, kk :== k), .])
        for (b=1; b<=B; b++) {
            hi = (b < B ? edges[b+1] : 1e9)
            S[k, b] = mean(Lk :>= edges[b] :& Lk :< hi)
        }
    }
    small = rowsum(S[., 1..3])
    // thresholds: distinct levels in random order, one change for the first nchange countries
    tau  = J(C, T, .)
    tau0 = J(C, 1, .)
    perm = order(runiform(C, 1), 1)
    for (c=1; c<=C; c++) {
        t0  = grid[mod(perm[c]-1, 5) + 1]
        tc  = 2 + ceil(4*runiform(1,1))
        t1  = t0 + (runiform(1,1) < .5 ? -.5 : .5)
        if (t1 >  1) t1 = t0 - .5
        if (t1 < -1) t1 = t0 + .5
        tau0[c] = t0
        for (t=1; t<=T; t++) {
            tau[c, t] = (c <= nchange & t >= tc) ? t1 : t0
        }
    }
    // the omitted factor: the share of firms in the size bin just above the initial threshold
    Wb = J(K, C, .)
    for (c=1; c<=C; c++) {
        bfirst = select(1::B, edges' :== tau0[c])
        for (k=1; k<=K; k++) Wb[k, c] = S[k, bfirst]
    }
    wbar = mean(vec(Wb))
    // cells
    CT = C*T; N = K*CT
    x = J(N,1,.); z = J(N,1,.); zlev = J(N,1,.); y = J(N,1,.); Zb = J(N,B,.)
    ctry = J(N,1,.); ind = J(N,1,.); cy = J(N,1,.); kc = J(N,1,.); iy = J(N,1,.)
    alpha = rnormal(K, 1, 0, .5)
    delta = rnormal(CT, 1, 0, .5)
    rho   = rnormal(C, 1, 0, sig_rho)
    for (k=1; k<=K; k++) {
        for (c=1; c<=C; c++) {
            i0  = (k-1)*C + c
            u   = rnormal(T, 1, 0, sig_u)
            th  = tau[c, .]'
            Lt  = J(T, 1, 1) * L0[i0, .] :+ lambda :* u :+ rnormal(T, n, 0, sig_nu)
            comply = (runiform(T, n) :< e)
            xs  = rowsum((Lt :>= th) :& comply) :/ n
            g0  = (edges :>= tau0[c])
            for (t=1; t<=T; t++) {
                i  = (k-1)*CT + (c-1)*T + t
                ct = (c-1)*T + t
                g  = (edges :>= tau[c, t])
                Zb[i, .] = S[k, .] :* g
                z[i]    = sum(Zb[i, .])
                zlev[i] = sum(S[k, .] :* g0)
                x[i]    = xs[t]
                y[i]    = alpha[k] + delta[ct] + beta*x[i] + u[t] ///
                        + gamma*(Wb[k, c] - wbar) + rho[c]*small[k] + rnormal(1,1,0,sig_e)
                ctry[i] = c; ind[i] = k; cy[i] = ct; kc[i] = i0; iy[i] = (k-1)*T + t
            }
        }
    }
}

// ---------- estimators on one draw ---------------------------------------------
// Returns (b_ols, b_rf, b_iv, F, pi, se_cy, se_rf_cy, se_2way, se_country):
// standard errors of 2SLS clustered by country-year, of the reduced form by
// country-year, of 2SLS two-way by country-industry and country-year, and of
// 2SLS by country. Fills the Rotemberg weights by bin and the country weights.
real rowvector est(real colvector x, real colvector z, real colvector y,
                   real matrix Zb, real colvector ctry, real scalar K,
                   real scalar C, real scalar T, real scalar mode, real colvector cy,
                   real colvector kc,
                   real colvector ab, real colvector bb,
                   real colvector ac, real colvector bc)
{
    N  = rows(x); CT = C*T; KC = K*C
    xp = dm(x, K, C, T, mode); yp = dm(y, K, C, T, mode); zp = dm(z, K, C, T, mode)
    bols = (xp'yp)/(xp'xp)
    brf  = (zp'yp)/(zp'zp)
    biv  = (zp'yp)/(zp'xp)
    pi   = (zp'xp)/(zp'zp)
    ex   = xp - pi*zp
    df   = (mode == 1 ? N - K*T - CT + T - 1 : N - (KC + K*T + CT - K - C - T + 1) - 1)
    se   = sqrt((ex'ex/df)/(zp'zp))
    F    = (pi/se)^2
    // clustered standard errors
    eiv = yp - biv*xp
    erf = yp - brf*zp
    giv = zp :* eiv
    grf = zp :* erf
    s_cy = 0; r_cy = 0; s_kc = 0; s_c = 0
    for (g=1; g<=CT; g++) {
        sel = (cy :== g)
        s_cy = s_cy + (sum(giv :* sel))^2
        r_cy = r_cy + (sum(grf :* sel))^2
    }
    for (g=1; g<=KC; g++) {
        sel = (kc :== g)
        s_kc = s_kc + (sum(giv :* sel))^2
    }
    for (g=1; g<=C; g++) {
        sel = (ctry :== g)
        s_c = s_c + (sum(giv :* sel))^2
    }
    s_hc  = sum(giv :^ 2)
    d2    = (zp'xp)^2
    se_cy = sqrt(CT/(CT-1)*s_cy/d2)
    se_rf = sqrt(CT/(CT-1)*r_cy)/(zp'zp)
    v2w   = KC/(KC-1)*s_kc + CT/(CT-1)*s_cy - s_hc
    se_2w = sqrt(max((v2w, 0))/d2)
    se_c  = sqrt(C/(C-1)*s_c/d2)
    // Rotemberg weights by size bin
    B  = cols(Zb)
    ab = J(B,1,0); bb = J(B,1,.)
    for (b=1; b<=B; b++) {
        zb = dm(Zb[., b], K, C, T, mode)
        d  = zb'xp
        ab[b] = d/(zp'xp)
        bb[b] = (abs(d) < 1e-12 ? . : (zb'yp)/d)
    }
    // the same decomposition by country
    ac = J(C,1,0); bc = J(C,1,.)
    for (c=1; c<=C; c++) {
        zc = dm(z :* (ctry :== c), K, C, T, mode)
        d  = zc'xp
        ac[c] = d/(zp'xp)
        bc[c] = (abs(d) < 1e-12 ? . : (zc'yp)/d)
    }
    return((bols, brf, biv, F, pi, se_cy, se_rf, se_2w, se_c))
}
end

*==============================================================================
* 1. One draw: the design, the identities, and the decompositions
*==============================================================================
mata:
K = 30; C = 10; T = 8; n = 50
x = z = zlev = y = Zb = ctry = ind = cy = kc = iy = S = small = tau = Wb = .
ab = bb = ac = bc = .
rseed(2025)
draw(K, C, T, n, 0.5, 0, 0, 3, 1, `beta', x, z, zlev, y, Zb, ctry, ind, cy, kc, iy, S, small, tau, Wb)
"One draw (lambda = 0.5, gamma = 0, sig_rho = 0, e = 1; 3 of 10 countries change their threshold)"
"  Structure 1 (industry-year and country-year effects, the paper):"
"  OLS, RF, 2SLS, F, first stage, SE 2SLS (cy), SE RF (cy), SE 2SLS (two-way), SE 2SLS (country)"
r1 = est(x, z, y, Zb, ctry, K, C, T, 1, cy, kc, ab, bb, ac, bc)
ab1 = ab; bb1 = bb; ac1 = ac; bc1 = bc
r1
"  Structure 2 (plus country-industry effects):"
r2 = est(x, z, y, Zb, ctry, K, C, T, 2, cy, kc, ab, bb, ac, bc)
r2
"  Rotemberg weights by size bin (bins cut at -1, -0.5, 0, 0.5, 1) and their estimates, structure 1:"
(ab1, bb1)'
"  Country weights and their estimates, structure 1:"
(ac1, bc1)'
"  Country weights and their estimates, structure 2 (only countries whose threshold changes carry weight):"
(ac, bc)'
"  sums of bin and country weights, structures 1 and 2:"
(sum(ab1), sum(ac1), sum(ab), sum(ac))
"  weighted bin estimates, weighted country estimates, 2SLS (structure 1, then 2):"
(sum(select(ab1, bb1 :< .) :* select(bb1, bb1 :< .)), sum(select(ac1, bc1 :< .) :* select(bc1, bc1 :< .)), r1[3])
(sum(select(ab, bb :< .) :* select(bb, bb :< .)), sum(select(ac, bc :< .) :* select(bc, bc :< .)), r2[3])
assert(abs(sum(ab1) - 1) < 1e-8)
assert(abs(sum(ac1) - 1) < 1e-8)
assert(abs(sum(select(ab1, bb1 :< .) :* select(bb1, bb1 :< .)) - r1[3]) < 1e-8)
assert(abs(sum(select(ac1, bc1 :< .) :* select(bc1, bc1 :< .)) - r1[3]) < 1e-8)
assert(abs(sum(ab) - 1) < 1e-8)
assert(abs(sum(ac) - 1) < 1e-8)
assert(abs(sum(select(ab, bb :< .) :* select(bb, bb :< .)) - r2[3]) < 1e-8)
assert(abs(sum(select(ac, bc :< .) :* select(bc, bc :< .)) - r2[3]) < 1e-8)
// which variation identifies: the threshold levels or the threshold changes
z1 = dm(z, K, C, T, 1); l1 = dm(zlev, K, C, T, 1); d1 = z1 - l1
l2 = dm(zlev, K, C, T, 2)
"  Structure 1: variance of demeaned z, of its level part, of its change part, and 2 Cov:"
(variance(z1), variance(l1), variance(d1), 2*(l1'd1)/(rows(z1)-1))
levR2 = ((l1'z1)/(rows(z1)-1))^2/(variance(l1)*variance(z1))
"  Structure 1: R-squared of demeaned z on its level part:"
levR2
"  Structure 2: max abs demeaned level part (should be zero):"
max(abs(l2))
assert(max(abs(l2)) < 1e-10)
// the omitted factor's correlation with the identifying variation
wcell = J(rows(z), 1, .)
for (k=1; k<=K; k++) {
    for (c=1; c<=C; c++) {
        for (t=1; t<=T; t++) wcell[(k-1)*C*T + (c-1)*T + t] = Wb[k, c]
    }
}
w1 = dm(wcell, K, C, T, 1)
"  Structure 1: correlation of the demeaned omitted factor with demeaned z:"
corrWZ = (w1'z1)/sqrt((w1'w1)*(z1'z1))
corrWZ
// base shares of three example industries for the deck's threshold stage
ex_k = (3, 15, 28)
"  base shares of industries 3, 15, 28 (small, middle, large firms):"
S[ex_k, .]
st_addobs(rows(x))
idx = st_addvar("double", ("x", "z", "y", "ind", "cy", "ctry", "kc", "iy"))
st_store(., idx, (x, z, y, ind, cy, ctry, kc, iy))
st_matrix("r1", r1); st_matrix("r2", r2)
st_numscalar("levR2", levR2); st_numscalar("corrWZ", corrWZ)
st_matrix("S_ex", S[ex_k, .]); st_matrix("ab", ab1); st_matrix("bb", bb1); st_matrix("ac", ac1); st_matrix("bc", bc1)
st_matrix("ac2", ac); st_matrix("bc2", bc)
st_matrix("tau", tau)
end

label var x "share of firms above the threshold that report (scope of the mandate)"
label var z "shift-share measure: base-period share above the threshold"
label var y "innovation in the cell"
di _n "Stata's own estimates on the same draw, structure 1 (industry-year and country-year effects)"
regress y x i.iy i.cy, vce(cluster cy)
scalar s1_ols = _b[x]
regress y z i.iy i.cy, vce(cluster cy)
scalar s1_rf = _b[z]
scalar s1_rfse = _se[z]
ivregress 2sls y i.iy i.cy (x = z), vce(cluster cy)
scalar s1_iv = _b[x]
scalar s1_se = _se[x]
estat firststage
quietly regress x z i.iy i.cy
scalar s1_fs = _b[z]
di _n "Structure 2 (plus country-industry effects)"
regress y x i.kc i.cy i.iy, vce(cluster cy)
scalar s2_ols = _b[x]
regress y z i.kc i.cy i.iy, vce(cluster cy)
scalar s2_rf = _b[z]
ivregress 2sls y i.kc i.cy i.iy (x = z), vce(cluster cy)
scalar s2_iv = _b[x]
scalar s2_se = _se[x]
estat firststage
quietly regress x z i.kc i.cy i.iy
scalar s2_fs = _b[z]
di _n "Identity checks (should be zero):"
di "  structure 1: regress OLS - demeaned x'y/x'x    = " %12.9f (s1_ols - r1[1,1])
di "  structure 1: regress RF  - demeaned z'y/z'z    = " %12.9f (s1_rf - r1[1,2])
di "  structure 1: ivregress   - demeaned z'y/z'x    = " %12.9f (s1_iv - r1[1,3])
di "  structure 1: 2SLS - reduced form / first stage  = " %12.9f (s1_iv - s1_rf/s1_fs)
di "  structure 2: regress OLS - demeaned x'y/x'x    = " %12.9f (s2_ols - r2[1,1])
di "  structure 2: regress RF  - demeaned z'y/z'z    = " %12.9f (s2_rf - r2[1,2])
di "  structure 2: ivregress   - demeaned z'y/z'x    = " %12.9f (s2_iv - r2[1,3])
di "  structure 2: 2SLS - reduced form / first stage  = " %12.9f (s2_iv - s2_rf/s2_fs)
di "  clustered SE of 2SLS, ivregress / Mata (equal up to the small-sample factor): " %8.4f (s1_se/r1[1,6]) " and " %8.4f (s2_se/r2[1,6])
assert abs(s1_ols - r1[1,1]) < 1e-6
assert abs(s1_rf  - r1[1,2]) < 1e-6
assert abs(s1_iv  - r1[1,3]) < 1e-6
assert abs(s1_iv - s1_rf/s1_fs) < 1e-8
assert abs(s2_ols - r2[1,1]) < 1e-6
assert abs(s2_rf  - r2[1,2]) < 1e-6
assert abs(s2_iv  - r2[1,3]) < 1e-6
assert abs(s2_iv - s2_rf/s2_fs) < 1e-8
assert s1_se/r1[1,6] > 0.9 & s1_se/r1[1,6] < 1.2
assert s2_se/r2[1,6] > 0.9 & s2_se/r2[1,6] < 1.2
di _n "One draw, structure 1: OLS " %7.3f s1_ols ", reduced form " %7.3f s1_rf " (SE " %6.3f s1_rfse "), 2SLS " %7.3f s1_iv " (SE " %6.3f s1_se "), first stage " %6.3f s1_fs ", F " %8.1f r1[1,4]
di    "One draw, structure 2: OLS " %7.3f s2_ols ", reduced form " %7.3f s2_rf ", 2SLS " %7.3f s2_iv " (SE " %6.3f s2_se "), first stage " %6.3f s2_fs ", F " %8.1f r2[1,4]
di    "Causal effect `beta'. R-squared of the structure-1 demeaned z on its threshold-level part: " %6.3f levR2 "; correlation of the omitted factor with demeaned z: " %6.3f corrWZ

*==============================================================================
* 2. Monte Carlo grids (structure 1, the paper's, with structure 2 alongside)
*==============================================================================
mata:
rseed(20252)
R = `R'
lam_grid = (0, 0.05, 0.1, 0.2, 0.3, 0.5)
gam_grid = (0, 0.1, 0.25, 0.5, 1, 2)
C_grid   = (3, 5, 10, 20)
e_grid   = (1, 0.8, 0.6, 0.4, 0.2)
// Columns: 1 grid value; structure 1: 2 mean OLS, 3 mean RF, 4 mean 2SLS, 5-6 2SLS
// quartiles, 7 mean F, 8 mean first-stage slope, 9 sd of 2SLS across draws, 10 median
// largest absolute country weight, 11 mean SE clustered by country-year, 12 rejection
// rate of the true beta with that SE, 19 mean two-way clustered SE, 20 its rejection
// rate, 21 mean country-clustered SE, 22 its rejection rate; structure 2: 13 mean
// 2SLS, 14 sd, 15 mean RF, 16 mean F, 17-18 2SLS quartiles.
real matrix mc(real scalar R, real rowvector grid, real scalar which, real scalar K, real scalar C0, real scalar T, real scalar n, real scalar beta)
{
    x = z = zlev = y = Zb = ctry = ind = cy = kc = iy = S = small = tau = Wb = .
    ab = bb = ac = bc = .
    out = J(cols(grid), 22, .)
    for (j=1; j<=cols(grid); j++) {
        o  = J(R, 9, .)
        o2 = J(R, 9, .)
        w  = J(R, 1, .)
        for (r=1; r<=R; r++) {
            lambda = 0.5; gamma = 0; sig_rho = 0; nchange = 3; e = 1; C = C0
            if (which == 1) lambda = grid[j]
            if (which == 2) gamma  = grid[j]
            if (which == 3) {
                C = grid[j]
                sig_rho = 0.3
                nchange = min((3, C))
            }
            if (which == 4) e = grid[j]
            draw(K, C, T, n, lambda, gamma, sig_rho, nchange, e, beta, x, z, zlev, y, Zb, ctry, ind, cy, kc, iy, S, small, tau, Wb)
            o[r, .]  = est(x, z, y, Zb, ctry, K, C, T, 1, cy, kc, ab, bb, ac, bc)
            w[r]     = max(abs(ac))
            o2[r, .] = est(x, z, y, Zb, ctry, K, C, T, 2, cy, kc, ab, bb, ac, bc)
        }
        s   = sort(o[., 3], 1)
        s2  = sort(o2[., 3], 1)
        ws  = sort(w, 1)
        rej   = mean(abs(o[., 3] :- beta) :/ o[., 6] :> 1.96)
        rej2w = mean(abs(o[., 3] :- beta) :/ o[., 8] :> 1.96)
        rejc  = mean(abs(o[., 3] :- beta) :/ o[., 9] :> 1.96)
        out[j, .] = (grid[j], mean(o[., 1]), mean(o[., 2]), mean(o[., 3]), s[ceil(.25*R)], s[ceil(.75*R)], ///
                     mean(o[., 4]), mean(o[., 5]), sqrt(variance(o[., 3])), ws[ceil(.5*R)], mean(o[., 6]), rej, ///
                     mean(o2[., 3]), sqrt(variance(o2[., 3])), mean(o2[., 2]), mean(o2[., 4]), s2[ceil(.25*R)], s2[ceil(.75*R)], ///
                     mean(o[., 8]), rej2w, mean(o[., 9]), rejc)
    }
    return(out)
}
"Columns: grid, OLS, RF, 2SLS, q1, q3, F, first stage, sd 2SLS, median largest country weight, SE cy, rej cy (structure 1);"
"  2SLS, sd, RF, F, q1, q3 (structure 2); SE two-way, rej two-way, SE country, rej country (structure 1)"
"Grid (a): lambda (reverse causality)"
A = mc(R, lam_grid, 1, K, C, T, n, `beta')
A
"Grid (b): gamma (an omitted country-industry factor: the share of firms in the bin just above the initial threshold), lambda = 0.5"
G = mc(R, gam_grid, 2, K, C, T, n, `beta')
G
"Grid (c): number of countries, with exposure-correlated country factors (sig_rho = 0.3), lambda = 0.5"
H = mc(R, C_grid, 3, K, C, T, n, `beta')
H
"Grid (d): e (compliance with the mandate), lambda = 0.5"
Q = mc(R, e_grid, 4, K, C, T, n, `beta')
Q
st_matrix("A", A); st_matrix("G", G); st_matrix("H", H); st_matrix("Q", Q)
end
assert abs(A[1,4] - `beta') < 0.06
assert abs(A[3,4] - `beta') < 0.06
assert abs(G[1,4] - `beta') < 0.06
di _n "Checks: mean 2SLS at lambda = 0, lambda = 0.1, gamma = 0 within 0.06 of `beta'"

*==============================================================================
* 3. Write the deck's data file
*==============================================================================
tempname js
file open `js' using "`outjs'", write replace
file write `js' "/* Written by scripts/bartik_shift_share.do (seed 2025). Do not edit by hand. */" _n
file write `js' "window.SIMBARTIK = {" _n
file write `js' "  seed: 2025, K: 30, C: 10, T: 8, n: 50, R: `R', beta: `beta', nchange: 3," _n
file write `js' "  edges: [-1, -0.5, 0, 0.5, 1]," _n
file write `js' "  one: {" _n
file write `js' "    paperFE: { ols: " %8.4f (s1_ols) ", rf: " %8.4f (s1_rf) ", rfSE: " %8.4f (s1_rfse) ", iv: " %8.4f (s1_iv) ", ivSE: " %8.4f (s1_se) ", F: " %8.1f (r1[1,4]) ", pi: " %8.4f (s1_fs) ", levelR2: " %6.3f (levR2) ", corrWZ: " %6.3f (corrWZ) " }," _n
file write `js' "    withinFE: { ols: " %8.4f (s2_ols) ", rf: " %8.4f (s2_rf) ", iv: " %8.4f (s2_iv) ", ivSE: " %8.4f (s2_se) ", F: " %8.1f (r2[1,4]) ", pi: " %8.4f (s2_fs) " }," _n
file write `js' "    binWeights: ["
forvalues b = 1/6 {
    file write `js' %8.4f (ab[`b',1])
    if `b' < 6 file write `js' ", "
}
file write `js' "]," _n "    binEstimates: ["
forvalues b = 1/6 {
    if bb[`b',1] == . file write `js' "null"
    else file write `js' %8.4f (bb[`b',1])
    if `b' < 6 file write `js' ", "
}
file write `js' "]," _n "    countryWeights: ["
forvalues c = 1/10 {
    file write `js' %8.4f (ac[`c',1])
    if `c' < 10 file write `js' ", "
}
file write `js' "]," _n "    countryEstimates: ["
forvalues c = 1/10 {
    if bc[`c',1] == . file write `js' "null"
    else file write `js' %8.4f (bc[`c',1])
    if `c' < 10 file write `js' ", "
}
file write `js' "]," _n "    countryWeightsWithin: ["
forvalues c = 1/10 {
    file write `js' %8.4f (ac2[`c',1])
    if `c' < 10 file write `js' ", "
}
file write `js' "]," _n "    exampleShares: ["
forvalues k = 1/3 {
    file write `js' "["
    forvalues b = 1/6 {
        file write `js' %7.4f (S_ex[`k',`b'])
        if `b' < 6 file write `js' ", "
    }
    file write `js' "]"
    if `k' < 3 file write `js' ", "
}
file write `js' "]," _n "    tau: ["
forvalues c = 1/10 {
    file write `js' "["
    forvalues t = 1/8 {
        file write `js' %5.2f (tau[`c',`t'])
        if `t' < 8 file write `js' ", "
    }
    file write `js' "]"
    if `c' < 10 file write `js' ", "
}
file write `js' "] }," _n
foreach g in A G H Q {
    if "`g'" == "A" {
        local key "lambda"
        file write `js' "  lambda: [" _n
    }
    if "`g'" == "G" {
        local key "gamma"
        file write `js' "  gamma: [" _n
    }
    if "`g'" == "H" {
        local key "C"
        file write `js' "  countries: [" _n
    }
    if "`g'" == "Q" {
        local key "e"
        file write `js' "  compliance: [" _n
    }
    local rows = rowsof(`g')
    forvalues j = 1/`rows' {
        file write `js' "    { `key': " %6.2f (`g'[`j',1]) ", ols: " %8.4f (`g'[`j',2]) ", rf: " %8.4f (`g'[`j',3]) ", iv: " %8.4f (`g'[`j',4]) ///
            ", q1: " %8.4f (`g'[`j',5]) ", q3: " %8.4f (`g'[`j',6]) ", F: " %9.1f (`g'[`j',7]) ", pi: " %8.4f (`g'[`j',8]) ", sd: " %8.4f (`g'[`j',9]) ", maxw: " %8.4f (`g'[`j',10]) ///
            ", seCY: " %8.4f (`g'[`j',11]) ", rejCY: " %6.3f (`g'[`j',12]) ", se2W: " %8.4f (`g'[`j',19]) ", rej2W: " %6.3f (`g'[`j',20]) ", seC: " %8.4f (`g'[`j',21]) ", rejC: " %6.3f (`g'[`j',22]) ///
            ", ivWithin: " %8.4f (`g'[`j',13]) ", sdWithin: " %8.4f (`g'[`j',14]) ", rfWithin: " %8.4f (`g'[`j',15]) ", FWithin: " %9.1f (`g'[`j',16]) ", q1Within: " %8.4f (`g'[`j',17]) ", q3Within: " %8.4f (`g'[`j',18]) " }"
        if `j' < `rows' file write `js' ","
        file write `js' _n
    }
    file write `js' "  ]," _n
}
file write `js' "};" _n
file close `js'
di _n "Done: `outjs'"

*==============================================================================
* 4. Figures (twoway, base Stata)
*==============================================================================
local cols "ols rf iv q1 q3 F pi sd maxw seCY rejCY ivW sdW rfW FW q1W q3W se2W rej2W seC rejC"
clear
svmat A, names(col)
rename (c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22) (lambda `cols')
twoway (line ols lambda, lcolor("144 53 59") lwidth(medthick)) ///
       (line rf lambda, lcolor("26 71 111") lwidth(medthick) lpattern(dash)) ///
       (line iv lambda, lcolor("26 71 111") lwidth(medthick)) ///
       (line q1 lambda, lcolor("26 71 111") lpattern(shortdash)) ///
       (line q3 lambda, lcolor("26 71 111") lpattern(shortdash)), ///
       yline(`beta', lcolor(gs8) lpattern(dash)) ///
       ytitle("estimate") xtitle("lambda (reverse causality)") ///
       title("OLS drifts with reverse causality; shift-share does not") ///
       legend(order(1 "OLS on the actual share" 2 "reduced form on the shift-share measure" 3 "2SLS" 4 "2SLS interquartile range") rows(2) position(6)) ///
       graphregion(color(white)) plotregion(margin(small))
graph export "figures/fig1_reverse_causality.png", replace width(2000)
clear
svmat G, names(col)
rename (c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22) (gamma `cols')
twoway (line iv gamma, lcolor("26 71 111") lwidth(medthick)) ///
       (line q1 gamma, lcolor("26 71 111") lpattern(shortdash)) ///
       (line q3 gamma, lcolor("26 71 111") lpattern(shortdash)) ///
       (line ivW gamma, lcolor(gs6) lwidth(medthick) lpattern(dash)) ///
       (line q1W gamma, lcolor(gs10) lpattern(shortdash)) ///
       (line q3W gamma, lcolor(gs10) lpattern(shortdash)), ///
       yline(`beta', lcolor(gs8) lpattern(dash)) ///
       ytitle("2SLS estimate") xtitle("gamma (omitted country-industry factor)") ///
       title("An omitted country-industry factor biases the design") ///
       legend(order(1 "country-year and industry-year effects (the paper)" 2 "its interquartile range" 4 "plus country-industry effects" 5 "its interquartile range") rows(2) position(6)) ///
       graphregion(color(white)) plotregion(margin(small))
graph export "figures/fig2_exclusion_violation.png", replace width(2000)
clear
svmat H, names(col)
rename (c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22) (C `cols')
twoway (connected sd C, lcolor("26 71 111") mcolor("26 71 111") lwidth(medthick)) ///
       (connected se2W C, lcolor("144 53 59") mcolor("144 53 59") lwidth(medthick) lpattern(dash)) ///
       (connected seC C, lcolor(gs6) mcolor(gs6) lwidth(medthick) lpattern(shortdash)) ///
       (connected maxw C, lcolor("184 134 11") mcolor("184 134 11") lwidth(medthick) yaxis(2)), ///
       ytitle("SD of 2SLS across draws; mean clustered SE") ytitle("largest country weight", axis(2)) ///
       xtitle("number of countries") xlabel(3 5 10 20) ///
       title("Few countries: one country dominates") ///
       legend(order(1 "SD of 2SLS across draws" 2 "mean SE, two-way clustered (the paper's scheme)" 3 "mean SE clustered by country" 4 "largest country weight (right axis)") rows(2) position(6)) ///
       graphregion(color(white)) plotregion(margin(small))
graph export "figures/fig3_few_countries.png", replace width(2000)
clear
svmat Q, names(col)
rename (c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22) (e `cols')
twoway (line pi e, lcolor(gs6) lwidth(medthick) lpattern(dot)) ///
       (line rf e, lcolor("26 71 111") lwidth(medthick) lpattern(dash)) ///
       (line iv e, lcolor("26 71 111") lwidth(medthick)), ///
       yline(`beta', lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs12)) ///
       ytitle("estimate") xtitle("e (compliance with the mandate)") ///
       title("Low compliance shrinks the reduced form, not 2SLS") ///
       legend(order(1 "first-stage slope of x on z" 2 "reduced form" 3 "2SLS") rows(1) position(6)) ///
       graphregion(color(white)) plotregion(margin(small))
graph export "figures/fig4_compliance.png", replace width(2000)
di _n "Done: figures/fig1_reverse_causality.png, fig2_exclusion_violation.png, fig3_few_countries.png, fig4_compliance.png"
