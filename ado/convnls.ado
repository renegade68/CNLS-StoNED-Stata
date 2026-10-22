*! version 0.5.2  07jul2026
*! convnls - Convex Nonparametric Least Squares (StoNED Stage 1)
*!
*! Kuosmanen (2008, 2012) CNLS with Cutting Plane + ADMM.
*!
*! SYNTAX: convnls depvar [if] [in], INputs(varlist) [options]
*!
*! Author: Choonjoo Lee (bloom.rampike@gmail.com)

program define convnls, eclass sortpreserve
    version 16.0
    syntax varlist(max=1 numeric) [if] [in], ///
        INputs(varlist numeric min=1) ///
        [FRONTier(string)] ///
        [GENerate(name)] ///
        [FITted(name)] ///
        [CONtext(varlist numeric)] ///
        [RTS(string)] ///
        [noMONotone] ///
        [RHO(real -1)] ///
        [TOLerance(real 1e-4)] ///
        [MAXiter(integer 20000)] ///
        [MAXCut(integer 100)] ///
        [ADDper(integer 0)] ///
        [PRUne(real 0)] ///
        [INit(string)] ///
        [ADAMIter(integer 100)] ///
        [ADAMLR(real 0.01)] ///
        [ADAMLambda(real 10)] ///
        [ADAMScreen(integer 2000)] ///
        [SWEEtpct(real 3)] ///
        [VCstrat(string)] ///
        [VCE(string)] ///
        [VERBose] ///
        [Level(cilevel)]

    marksample touse
    local depvar `varlist'

    quietly count if `touse'
    if r(N) == 0 error 2000
    local n_obs = r(N)

    local m : word count `inputs'
    local q : word count `context'

    if "`frontier'" == "" local frontier "production"
    local frontier = lower("`frontier'")
    if !inlist("`frontier'", "production", "cost") {
        di as error "frontier() must be {bf:production} or {bf:cost}"
        exit 198
    }
    local is_concave = ("`frontier'" == "production")

    if "`rts'" == "" local rts "vrs"
    local rts = lower("`rts'")
    if !inlist("`rts'", "crs", "vrs") {
        di as error "rts() must be {bf:crs} or {bf:vrs}"
        exit 198
    }
    if `q' > 0 & "`rts'" == "crs" {
        di as error "context() is currently supported with rts(vrs) only"
        exit 198
    }

    if "`vce'" == "" local vce "ols"
    local vce = lower("`vce'")
    if !inlist("`vce'", "ols", "robust") {
        di as error "vce() must be {bf:ols} or {bf:robust}"
        exit 198
    }

    local do_mono = 1
    if "`monotone'" == "nomonotone" local do_mono = 0

    if "`init'" == "" local init "none"
    local init = lower("`init'")
    if !inlist("`init'", "none", "ele", "sweet", "adam") {
        di as error "init() must be {bf:none}, {bf:ele}, {bf:sweet}, or {bf:adam}"
        exit 198
    }
    if `adamiter' < 0 {
        di as error "adamiter() must be nonnegative"
        exit 198
    }
    if `adamlr' <= 0 {
        di as error "adamlr() must be positive"
        exit 198
    }
    if `adamlambda' <= 0 {
        di as error "adamlambda() must be positive"
        exit 198
    }
    if `adamscreen' < 0 {
        di as error "adamscreen() must be nonnegative"
        exit 198
    }

    if "`vcstrat'" == "" local vcstrat "all"
    local vcstrat = lower("`vcstrat'")
    if !inlist("`vcstrat'", "all", "group", "one") {
        di as error "vcstrat() must be {bf:all}, {bf:group}, or {bf:one}"
        exit 198
    }

    // Avoid premature termination when vcstrat(one) uses the default maxcut.
    if "`vcstrat'" == "one" & `maxcut' == 100 {
        local maxcut = min(`n_obs' * (`n_obs' - 1), 5000)
        di as text "Note: vcstrat(one) — maxcut auto-set to " `maxcut'
    }

    if `addper' < 0 {
        di as error "addper() must be nonnegative"
        exit 198
    }

    if "`generate'" == "" local generate "_convnls_resid"
    if "`fitted'" == ""   local fitted "_convnls_fitted"

    foreach v in `generate' `fitted' _convnls_alpha {
        capture drop `v'
    }
    forvalues j = 1/`m' {
        capture drop _convnls_beta`j'
    }

    foreach v of varlist `depvar' `inputs' `context' {
        quietly count if missing(`v') & `touse'
        if r(N) > 0 {
            di as error "Variable `v' has `r(N)' missing values"
            exit 459
        }
    }

    quietly gen double `generate' = .
    quietly gen double `fitted' = .
    quietly gen double _convnls_alpha = .
    forvalues j = 1/`m' {
        quietly gen double _convnls_beta`j' = .
    }

    di as text ""
    di as text "{hline 65}"
    di as text "CNLS: Convex Nonparametric Least Squares (StoNED Stage 1)"
    di as text "{hline 65}"
    di as text "Observations:      " as result %6.0f `n_obs'
    di as text "Output variable:   " as result "`depvar'"
    di as text "Input variables:   " as result "`inputs' (`m')"
    if `q' > 0 {
        di as text "Context variables: " as result "`context' (`q')"
        di as text "Context VCE:       " as result upper("`vce'")
    }
    di as text "Frontier type:     " as result ///
        cond(`is_concave', "Production (concave)", "Cost (convex)")
    di as text "Returns to scale:  " as result upper("`rts'")
    di as text "Monotonicity:      " as result cond(`do_mono', "Yes", "No")
    di as text "Tolerance:         " as result %9.1e `tolerance'
    if `addper' > 0 {
        di as text "Add per iter:      " as result %9.0f `addper' ///
            as text " largest violations"
    }
    if "`init'" == "adam" {
        di as text "Adam warm start:   " as result `adamiter' ///
            as text " iter, screen " as result `adamscreen'
    }
    di as text "{hline 65}"

    mata: _convnls_estimate(       ///
        "`depvar'", "`inputs'", "`context'", ///
        "`touse'",              ///
        "`generate'", "`fitted'", ///
        "`rts'", `do_mono', `is_concave', ///
        `rho', `tolerance',     ///
        `maxiter', `maxcut', `addper', `prune', ///
        "`init'", `sweetpct', "`vcstrat'", "`vce'", ///
        `adamiter', `adamlr', `adamlambda', `adamscreen', ///
        "`verbose'")

    // Save Mata r() results before any r-class command overwrites them.
    tempname mata_cp_history
    local mata_converged = 0
    local mata_cp_iters  = .
    local mata_n_active  = .
    local mata_max_viol  = .
    capture {
        local mata_converged = r(converged)
        local mata_cp_iters  = r(cp_iters)
        local mata_n_active  = r(n_active)
        local mata_max_viol  = r(max_viol)
        matrix `mata_cp_history' = r(cp_history)
    }
    tempname context_b context_V context_V_robust
    if `q' > 0 {
        matrix `context_b' = r(context_b)
        matrix `context_V' = r(context_V)
        matrix `context_V_robust' = r(context_V_robust)
        matrix colnames `context_b' = `context'
        matrix colnames `context_V' = `context'
        matrix rownames `context_V' = `context'
        matrix colnames `context_V_robust' = `context'
        matrix rownames `context_V_robust' = `context'
        if "`vce'" == "robust" {
            matrix `context_V' = `context_V_robust'
            matrix colnames `context_V' = `context'
            matrix rownames `context_V' = `context'
        }
    }

    tempvar resid_sq
    quietly gen double `resid_sq' = `generate'^2 if `touse'
    quietly summ `resid_sq' if `touse', meanonly
    local rss = r(sum)
    local sigma2 = `rss' / `n_obs'

    di as text ""
    di as text "Results:"
    di as text "{hline 40}"
    di as text "RSS:               " as result %12.4f `rss'
    di as text "sigma^2 (RSS/N):   " as result %12.6f `sigma2'
    di as text "sigma   (RMSE):    " as result %12.6f sqrt(`sigma2')

    quietly summ `generate' if `touse', detail
    local skew = r(skewness)
    di as text "Residual mean:     " as result %12.6f r(mean)
    di as text "Residual skewness: " as result %12.4f `skew'
    di as text "Residual kurtosis: " as result %12.4f r(kurtosis)

    quietly summ `fitted' if `touse', meanonly
    di as text "Fitted range:      [" as result %8.2f r(min) ///
        as text ", " as result %8.2f r(max) as text "]"

    if `is_concave' {
        if `skew' > 0 {
            di as text ""
            di as text "{bf:Note:} Residual skewness > 0 (" ///
                as result %5.3f `skew' as text ")."
            di as text "  Production frontier expects negative skewness."
        }
        else {
            di as text ""
            di as text "Residual skewness < 0 (" as result %5.3f `skew' ///
                as text ") -- consistent with v - u."
        }
    }
    else {
        if `skew' < 0 {
            di as text ""
            di as text "{bf:Note:} Residual skewness < 0."
            di as text "  Cost frontier expects positive skewness."
        }
        else {
            di as text ""
            di as text "Residual skewness > 0 (" as result %5.3f `skew' ///
                as text ") -- consistent with v + u."
        }
    }
    di as text "{hline 65}"

    label variable `generate' "CNLS residual"
    label variable `fitted' "CNLS fitted value"
    label variable _convnls_alpha "CNLS intercept (alpha_i)"
    forvalues j = 1/`m' {
        local vname : word `j' of `inputs'
        label variable _convnls_beta`j' "CNLS slope for `vname'"
    }

    ereturn clear
    if `q' > 0 {
        ereturn post `context_b' `context_V', esample(`touse')
    }
    ereturn scalar N = `n_obs'
    ereturn scalar m = `m'
    ereturn scalar q = `q'
    ereturn scalar rss = `rss'
    ereturn scalar sigma2 = `sigma2'
    ereturn scalar converged = `mata_converged'
    ereturn scalar cp_iters  = `mata_cp_iters'
    ereturn scalar n_active  = `mata_n_active'
    ereturn scalar max_viol  = `mata_max_viol'
    ereturn scalar addper    = `addper'
    ereturn scalar adamiter  = `adamiter'
    ereturn scalar adamlr    = `adamlr'
    ereturn scalar adamlambda = `adamlambda'
    ereturn scalar adamscreen = `adamscreen'
    capture ereturn matrix cp_history = `mata_cp_history'

    ereturn local depvar     "`depvar'"
    ereturn local inputs     "`inputs'"
    ereturn local context    "`context'"
    ereturn local frontier   "`frontier'"
    ereturn local rts        "`rts'"
    ereturn local resid_var    "`generate'"
    ereturn local fitted_var   "`fitted'"
    ereturn local init_method  "`init'"
    ereturn local vcstrat      "`vcstrat'"
    ereturn local vce          "`vce'"
    ereturn local cmd          "convnls"
    if `q' > 0 {
        ereturn matrix V_robust = `context_V_robust'
        di as text ""
        di as text "Context-variable coefficients (one-stage CNLS-z):"
        ereturn display, level(`level')
    }
end


// ============================================================================
// MATA ENGINE
// ============================================================================
version 16.0
mata:
mata set matastrict on

// ============================================================================
// Main estimation
// ============================================================================
void _convnls_estimate(string scalar depvar_str, string scalar inputs_str,
                    string scalar context_str,
                    string scalar touse_str,
                    string scalar resid_var, string scalar fitted_var,
                    string scalar rts, real scalar do_mono,
                    real scalar is_concave,
                    real scalar rho_init, real scalar tol,
                    real scalar maxiter, real scalar maxcut,
                    real scalar addper, real scalar prune_tol,
                    string scalar init_str, real scalar sweetpct_val,
                    string scalar vcstrat_str, string scalar vce_str,
                    real scalar adamiter, real scalar adamlr,
                    real scalar adamlambda, real scalar adamscreen,
                    string scalar verbose_str)
{
    real matrix    X, X_raw, C, C_raw, C_center, ZtZ, G_mono, G_afriat, G_full
    real matrix    beta_curr, violations, new_rows, cp_history
    real colvector y, y_raw, obs_idx, Zty, delta, alpha_curr
    real colvector resid, resid_scaled, fitted, x_scale, c_mean, c_scale
    real colvector theta_scaled, theta_orig
    real matrix    context_V, context_V_robust, XtX_inv, meat
    real rowvector zi, ci
    real scalar    N, m, q, p, p_local, verbose, i, j, base, cbase
    real scalar    cp_iter, max_viol, n_active_afriat, row_mono
    real scalar    resid_idx, fitted_idx, alpha_idx, beta_j_idx
    real scalar    n_mono_viol, rho, inner_tol
    real scalar    y_mean, y_sd, n_pruned, n_new_unique
    real scalar    near_converge_count, max_viol_hold, n_viol_hold
    // Constraint tracking: N x N matrix, 1 if (i,h) already added
    real matrix    added_pairs
    real scalar    is_crs
    real matrix    G_pairs
    real matrix    init_pairs
    real colvector delta_hold

    verbose = (verbose_str != "")

    // ---- Load data ----
    y_raw = st_data(., depvar_str, touse_str)
    X_raw = st_data(., inputs_str, touse_str)
    if (context_str == "") {
        C_raw = J(rows(X_raw), 0, .)
    }
    else {
        C_raw = st_data(., context_str, touse_str)
    }
    obs_idx = selectindex(st_data(., touse_str))

    N = rows(X_raw)
    m = cols(X_raw)
    q = cols(C_raw)
    is_crs = (rts == "crs")
    p_local = is_crs ? m : (1 + m)
    p = N * p_local + q

    if (verbose) {
        printf("{txt}Data: N=%g, m=%g, q=%g, local p=%g, total params=%g\n",
               N, m, q, p_local, p)
        printf("{txt}Potential Afriat constraints: %g\n", N*(N-1))
    }

    // ---- Data scaling ----
    y_mean = mean(y_raw)
    y_sd = sqrt(variance(y_raw))
    if (y_sd < 1e-10) y_sd = 1
    // CRS: no mean subtraction (alpha=0, so centering would add spurious intercept)
    if (is_crs) {
        y = y_raw :/ y_sd
    }
    else {
        y = (y_raw :- y_mean) :/ y_sd
    }

    x_scale = J(m, 1, 0)
    X = J(N, m, 0)
    for (j = 1; j <= m; j++) {
        x_scale[j] = sqrt(variance(X_raw[., j]))
        if (x_scale[j] < 1e-10) x_scale[j] = 1
        X[., j] = X_raw[., j] / x_scale[j]
    }
    C = J(N, q, 0)
    c_mean = J(q, 1, 0)
    c_scale = J(q, 1, 1)
    if (q > 0) {
        for (j = 1; j <= q; j++) {
            c_mean[j] = mean(C_raw[., j])
            c_scale[j] = sqrt(variance(C_raw[., j]))
            if (c_scale[j] < 1e-10) c_scale[j] = 1
            C[., j] = (C_raw[., j] :- c_mean[j]) / c_scale[j]
        }
    }

    if (verbose) {
        printf("{txt}Scaling: y_mean=%9.4f, y_sd=%9.4f\n", y_mean, y_sd)
    }

    // ---- Build Z'Z and Z'y ----
    ZtZ = J(p, p, 0)
    Zty = J(p, 1, 0)
    cbase = N * p_local

    for (i = 1; i <= N; i++) {
        if (is_crs) {
            zi = X[i, .]
        }
        else {
            zi = (1, X[i, .])
        }
        base = (i-1)*p_local
        ZtZ[|base+1, base+1 \ base+p_local, base+p_local|] =
            ZtZ[|base+1, base+1 \ base+p_local, base+p_local|] + zi' * zi
        Zty[|base+1 \ base+p_local|] = Zty[|base+1 \ base+p_local|] + y[i] * zi'
        if (q > 0) {
            ci = C[i, .]
            ZtZ[|base+1, cbase+1 \ base+p_local, cbase+q|] =
                ZtZ[|base+1, cbase+1 \ base+p_local, cbase+q|] + zi' * ci
            ZtZ[|cbase+1, base+1 \ cbase+q, base+p_local|] =
                ZtZ[|cbase+1, base+1 \ cbase+q, base+p_local|] + ci' * zi
            ZtZ[|cbase+1, cbase+1 \ cbase+q, cbase+q|] =
                ZtZ[|cbase+1, cbase+1 \ cbase+q, cbase+q|] + ci' * ci
            Zty[|cbase+1 \ cbase+q|] = Zty[|cbase+1 \ cbase+q|] + y[i] * ci'
        }
    }

    // ---- Rho initialization ----
    if (rho_init < 0) {
        rho = 1.0
    }
    else {
        rho = rho_init
    }

    // ---- OLS initialization ----
    delta = _convnls_init_ols(y, X, C, N, m, q, p_local, is_crs)
    if (verbose) printf("{txt}OLS init done, rho=%6.2f\n", rho)

    // ---- Monotonicity constraints ----
    if (do_mono) {
        G_mono = J(N*m, p, 0)
        row_mono = 0
        for (i = 1; i <= N; i++) {
            base = (i-1)*p_local
            for (j = 1; j <= m; j++) {
                row_mono++
                if (is_crs) {
                    G_mono[row_mono, base + j] = -1
                }
                else {
                    G_mono[row_mono, base + 1 + j] = -1
                }
            }
        }
    }
    else {
        G_mono = J(0, p, 0)
    }

    // ---- Initialize constraint tracking ----
    added_pairs = J(N, N, 0)

    // ---- Cutting Plane Main Loop ----
    G_afriat = J(0, p, 0)
    G_pairs  = J(0, 2, 0)
    n_active_afriat = 0
    near_converge_count = 0

    // ---- Init: pre-populate Afriat constraints (ele, sweet, or adam) ----
    if (init_str != "none") {
        if (verbose) printf("{txt}Init (%s): computing initial constraint set...\n", init_str)
        if (init_str == "ele") {
            init_pairs = _convnls_init_ele(X, N, added_pairs)
        }
        else if (init_str == "sweet") {
            init_pairs = _convnls_init_sweet(X, N, sweetpct_val, added_pairs)
        }
        else {   // adam
            delta = _convnls_init_adam(delta, y, X, ZtZ, Zty, N, m, p_local,
                                     is_crs, do_mono, is_concave,
                                     adamiter, adamlr, adamlambda,
                                     adamscreen, verbose)
            _convnls_extract(delta, N, m, p_local, is_crs, alpha_curr, beta_curr)
            init_pairs = _convnls_find_violations(alpha_curr, beta_curr, X,
                                                N, m, 0, is_concave)
            if (rows(init_pairs) > 0) {
                init_pairs = _convnls_sort_desc(init_pairs, 3)
                if (adamscreen > 0 & rows(init_pairs) > adamscreen) {
                    init_pairs = init_pairs[1..adamscreen, .]
                }
                _convnls_mark_added(init_pairs, added_pairs)
            }
        }
        if (rows(init_pairs) > 0) {
            new_rows        = _convnls_build_afriat_rows(init_pairs, X, N, m, p_local, q, is_concave, is_crs)
            G_afriat        = new_rows
            G_pairs         = init_pairs[., 1..2]
            n_active_afriat = rows(G_afriat)
            G_full = G_afriat \ G_mono
            if (verbose) {
                printf("{txt}  %g initial constraints; solving ADMM...\n", n_active_afriat)
            }
            delta = _convnls_admm_solve(ZtZ, Zty, G_full, delta, rho, tol, maxiter, verbose)
        }
    }

    printf("{txt}\nCutting Plane:\n")
    printf("{txt}  %4s  %7s  %7s  %7s  %10s\n",
           "Iter", "Active", "Violat", "NewAdd", "Max Viol")
    printf("{txt}  {hline 45}\n")
    cp_history = J(0, 5, .)

    for (cp_iter = 1; cp_iter <= maxcut; cp_iter++) {

        _convnls_extract(delta, N, m, p_local, is_crs, alpha_curr, beta_curr)

        violations = _convnls_find_violations(alpha_curr, beta_curr, X,
                                            N, m, tol, is_concave)

        if (rows(violations) == 0) {
            max_viol = 0
        }
        else {
            max_viol = max(violations[., 3])
        }

        n_pruned = 0
        if (prune_tol > 0 & n_active_afriat > 0 & cp_iter > 3) {
            _convnls_prune(G_afriat, G_pairs, added_pairs, delta, prune_tol,
                         N, m, p, n_pruned)
            n_active_afriat = rows(G_afriat)
        }

        if (max_viol <= tol) {
            printf("{txt}  %4.0f  %7.0f  %7.0f  %7s  %10.2e  *\n",
                   cp_iter, n_active_afriat, rows(violations), "--", max_viol)
            cp_history = cp_history \ (cp_iter, n_active_afriat, rows(violations), 0, max_viol)
            printf("{txt}  Converged.\n")
            st_numscalar("r(converged)", 1)
            st_numscalar("r(cp_iters)", cp_iter)
            st_numscalar("r(n_active)", n_active_afriat)
            st_numscalar("r(max_viol)", max_viol)
            st_matrix("r(cp_history)", cp_history)
            break
        }

        if (max_viol < 30 * tol) {
            near_converge_count = near_converge_count + 1
        }
        else {
            near_converge_count = 0
        }

        if (near_converge_count >= 2) {
            delta_hold    = delta
            max_viol_hold = max_viol
            n_viol_hold   = rows(violations)
            rho = max((rho, 1.0))
            G_full = G_afriat \ G_mono
            if (verbose) {
                printf("{txt}  Near-converged, final tight solve...\n")
            }
            delta = _convnls_admm_solve(ZtZ, Zty, G_full, delta,
                                      rho, tol * 0.1, maxiter, verbose)

            _convnls_extract(delta, N, m, p_local, is_crs, alpha_curr, beta_curr)
            violations = _convnls_find_violations(alpha_curr, beta_curr, X,
                                                N, m, tol, is_concave)
            if (rows(violations) == 0) {
                max_viol = 0
            }
            else {
                max_viol = max(violations[., 3])
            }

            printf("{txt}  %4.0f  %7.0f  %7.0f  %7s  %10.2e  **\n",
                   cp_iter, n_active_afriat, rows(violations), "final", max_viol)

            if (max_viol > max_viol_hold) {
                delta    = delta_hold
                max_viol = max_viol_hold
                printf("{txt}  Final solve worsened feasibility; retaining previous iterate (max_viol=%9.2e).\n", max_viol)
                cp_history = cp_history \ (cp_iter, n_active_afriat, n_viol_hold, -2, max_viol)
            }
            else {
                cp_history = cp_history \ (cp_iter, n_active_afriat, rows(violations), -1, max_viol)
            }

            if (max_viol <= tol) {
                printf("{txt}  Converged (after final solve).\n")
                st_numscalar("r(converged)", 1)
            }
            else {
                printf("{txt}  Practically converged (max_viol=%9.2e).\n", max_viol)
                st_numscalar("r(converged)", 2)
            }
            st_numscalar("r(cp_iters)", cp_iter)
            st_numscalar("r(n_active)", n_active_afriat)
            st_numscalar("r(max_viol)", max_viol)
            st_matrix("r(cp_history)", cp_history)
            break
        }

        if (vcstrat_str == "group") {
            violations = _convnls_filter_group(violations, N)
        }
        else if (vcstrat_str == "one") {
            if (rows(violations) > 0) {
                violations = _convnls_sort_desc(violations, 3)
                violations = violations[1, .]
            }
        }
        else {
            if (addper > 0 & rows(violations) > addper) {
                violations = _convnls_sort_desc(violations, 3)
                violations = violations[1..addper, .]
            }
        }

        n_new_unique = rows(violations)  // post-vcstrat, pre-dedup

        violations = _convnls_filter_new(violations, added_pairs)

        printf("{txt}  %4.0f  %7.0f  %7.0f  %7.0f  %10.2e\n",
               cp_iter, n_active_afriat, n_new_unique, rows(violations), max_viol)
        cp_history = cp_history \ (cp_iter, n_active_afriat, n_new_unique, rows(violations), max_viol)

        if (rows(violations) == 0) {
            // All violated constraints are already active; tighten feasibility.
            if (verbose) {
                printf("{txt}  No new constraints; boosting rho and re-solving\n")
            }
            rho = max((rho, 1.0))
            G_full = G_afriat \ G_mono
            delta = _convnls_admm_solve(ZtZ, Zty, G_full, delta,
                                      rho, tol * 0.1, maxiter, verbose)
            near_converge_count = near_converge_count + 1
            continue
        }

        // Mark pairs as added
        _convnls_mark_added(violations, added_pairs)

        // Build and append
        new_rows = _convnls_build_afriat_rows(violations, X, N, m, p_local, q, is_concave, is_crs)
        G_afriat = G_afriat \ new_rows
        G_pairs  = G_pairs  \ violations[., 1..2]
        n_active_afriat = rows(G_afriat)

        G_full = G_afriat \ G_mono

        // Adaptive inner tolerance
        inner_tol = max((0.01 * max_viol, tol))

        delta = _convnls_admm_solve(ZtZ, Zty, G_full, delta,
                                  rho, inner_tol, maxiter, verbose)

        if (cp_iter == maxcut) {
            printf("{txt}  Warning: max cutting plane iters\n")
            st_numscalar("r(converged)", 0)
            st_numscalar("r(cp_iters)", cp_iter)
            st_numscalar("r(n_active)", n_active_afriat)
            st_numscalar("r(max_viol)", max_viol)
            st_matrix("r(cp_history)", cp_history)
        }
    }

    // ---- Unscale and store ----
    _convnls_extract(delta, N, m, p_local, is_crs, alpha_curr, beta_curr)
    theta_scaled = J(q, 1, 0)
    theta_orig = J(q, 1, 0)
    if (q > 0) {
        theta_scaled = delta[|cbase+1 \ cbase+q|]
        for (j = 1; j <= q; j++) {
            theta_orig[j] = y_sd * theta_scaled[j] / c_scale[j]
        }
    }

    for (i = 1; i <= N; i++) {
        if (is_crs) {
            alpha_curr[i] = 0
        }
        else {
            alpha_curr[i] = y_mean + y_sd * alpha_curr[i]
            if (q > 0) {
                alpha_curr[i] = alpha_curr[i] - c_mean' * theta_orig
            }
        }
        for (j = 1; j <= m; j++) {
            beta_curr[i, j] = y_sd * beta_curr[i, j] / x_scale[j]
        }
    }

    fitted = J(N, 1, 0)
    for (i = 1; i <= N; i++) {
        fitted[i] = alpha_curr[i] + X_raw[i, .] * beta_curr[i, .]'
        if (q > 0) {
            fitted[i] = fitted[i] + C_raw[i, .] * theta_orig
        }
    }
    resid = y_raw - fitted
    resid_scaled = resid :/ y_sd

    if (q > 0) {
        C_center = C_raw
        for (j = 1; j <= q; j++) {
            C_center[., j] = C_raw[., j] :- c_mean[j]
        }
        XtX_inv = invsym(cross(C_center, C_center))
        context_V = (sum(resid:^2) / max((N - q, 1))) * XtX_inv
        meat = J(q, q, 0)
        for (i = 1; i <= N; i++) {
            meat = meat + (resid[i]^2) * C_center[i, .]' * C_center[i, .]
        }
        context_V_robust = (N / max((N - q, 1))) * XtX_inv * meat * XtX_inv
        st_matrix("r(context_b)", theta_orig')
        st_matrix("r(context_V)", context_V)
        st_matrix("r(context_V_robust)", context_V_robust)
    }

    resid_idx  = st_varindex(resid_var)
    fitted_idx = st_varindex(fitted_var)
    alpha_idx  = st_varindex("_convnls_alpha")

    for (i = 1; i <= N; i++) {
        st_store(obs_idx[i], resid_idx,  resid[i])
        st_store(obs_idx[i], fitted_idx, fitted[i])
        st_store(obs_idx[i], alpha_idx,  alpha_curr[i])
    }

    for (j = 1; j <= m; j++) {
        beta_j_idx = st_varindex("_convnls_beta" + strofreal(j))
        for (i = 1; i <= N; i++) {
            st_store(obs_idx[i], beta_j_idx, beta_curr[i, j])
        }
    }

    printf("{txt}\nFinal: %g active / %g possible (%5.1f%%)\n",
           n_active_afriat, N*(N-1), 100*n_active_afriat/(N*(N-1)))

    if (is_crs & n_active_afriat == 0) {
        printf("{txt}  Note: CRS fit ended with zero active Afriat constraints.\n")
        printf("{txt}  This is consistent with a global homogeneous linear fit;\n")
        printf("{txt}  check model fit if a curved CRS frontier was expected.\n")
    }

    if (do_mono) {
        n_mono_viol = 0
        for (i = 1; i <= N; i++) {
            for (j = 1; j <= m; j++) {
                if (beta_curr[i, j] < -tol * y_sd / x_scale[j]) n_mono_viol++
            }
        }
        if (n_mono_viol > 0) {
            printf("{txt}  Warning: %g monotonicity violations\n", n_mono_viol)
        }
        else {
            printf("{txt}  Monotonicity: OK\n")
        }
    }
}


// ============================================================================
// OLS initialization
// ============================================================================
real colvector _convnls_init_ols(real colvector y, real matrix X, real matrix C,
                               real scalar N, real scalar m, real scalar q,
                               real scalar p,
                               real scalar is_crs)
{
    real colvector delta, beta_ols, resid_ols, coef, xi
    real matrix    X_const, XX, Xy
    real scalar    alpha_ols, i, j, base, denom, scale_i, cbase

    if (is_crs) {
        delta = J(N*p + q, 1, 0)
        cbase = N*p
        for (i = 1; i <= N; i++) {
            base = (i-1)*p
            xi = X[i,.]'
            denom = cross(xi, xi)
            if (denom < 1e-10) denom = 1e-10
            scale_i = y[i] / denom
            for (j = 1; j <= m; j++) {
                delta[base + j] = max((scale_i * xi[j], 0))
            }
        }
    }
    else {
        X_const = J(N, 1, 1), X, C
        XX = cross(X_const, X_const)
        Xy = cross(X_const, y)
        coef = invsym(XX) * Xy

        alpha_ols = coef[1]
        beta_ols = coef[2..1+m]
        resid_ols = y - X_const * coef

        delta = J(N*p + q, 1, 0)
        cbase = N*p
        for (i = 1; i <= N; i++) {
            base = (i-1)*p
            delta[base + 1] = alpha_ols + resid_ols[i]
            for (j = 1; j <= m; j++) {
                delta[base + 1 + j] = max((beta_ols[j], 0))
            }
        }
        if (q > 0) {
            delta[|cbase+1 \ cbase+q|] = coef[|2+m \ 1+m+q|]
        }
    }
    return(delta)
}


// ============================================================================
// Adam warm start
// Penalized objective is used only to seed delta and active constraints.
// The final estimator remains the deterministic constrained CNLS solve.
// ============================================================================
real colvector _convnls_init_adam(real colvector delta0, real colvector y,
                               real matrix X, real matrix ZtZ,
                               real colvector Zty, real scalar N,
                               real scalar m, real scalar p,
                               real scalar is_crs, real scalar do_mono,
                               real scalar is_concave, real scalar maxiter,
                               real scalar lr, real scalar lambda,
                               real scalar screen, real scalar verbose)
{
    real colvector delta, grad, mom, vel
    real colvector alpha
    real matrix beta
    real matrix candidates
    real scalar t, i, h, j, k, base_i, base_h, idx_i, idx_h
    real scalar raw, viol, coeff, b1, b2, eps, step, bval, mv

    delta = delta0
    if (maxiter <= 0) return(delta)

    b1 = 0.9
    b2 = 0.999
    eps = 1e-8
    mom = J(rows(delta), 1, 0)
    vel = J(rows(delta), 1, 0)

    _convnls_extract(delta, N, m, p, is_crs, alpha, beta)
    candidates = _convnls_find_violations(alpha, beta, X, N, m, 0, is_concave)
    if (rows(candidates) > 0) {
        candidates = _convnls_sort_desc(candidates, 3)
        if (screen > 0 & rows(candidates) > screen) {
            candidates = candidates[1..screen, .]
        }
    }
    if (verbose) {
        printf("{txt}  Adam warm start candidates: %g\n", rows(candidates))
    }
    if (rows(candidates) == 0) return(delta)

    for (t = 1; t <= maxiter; t++) {
        grad = 2 * (ZtZ * delta - Zty)
        _convnls_extract(delta, N, m, p, is_crs, alpha, beta)

        for (k = 1; k <= rows(candidates); k++) {
            i = candidates[k, 1]
            h = candidates[k, 2]
            base_i = (i - 1) * p
            base_h = (h - 1) * p

            raw = alpha[i] - alpha[h] + (beta[i,.] - beta[h,.]) * X[h,.]'

            if (is_concave) {
                if (raw >= 0) continue
                viol = -raw
                coeff = -2 * lambda * viol
            }
            else {
                if (raw <= 0) continue
                viol = raw
                coeff = 2 * lambda * viol
            }

            if (!is_crs) {
                grad[base_i + 1] = grad[base_i + 1] + coeff
                grad[base_h + 1] = grad[base_h + 1] - coeff
            }
            for (j = 1; j <= m; j++) {
                idx_i = base_i + (is_crs ? j : 1 + j)
                idx_h = base_h + (is_crs ? j : 1 + j)
                grad[idx_i] = grad[idx_i] + coeff * X[h, j]
                grad[idx_h] = grad[idx_h] - coeff * X[h, j]
            }
        }

        if (do_mono) {
            for (i = 1; i <= N; i++) {
                base_i = (i - 1) * p
                for (j = 1; j <= m; j++) {
                    idx_i = base_i + (is_crs ? j : 1 + j)
                    bval = delta[idx_i]
                    if (bval < 0) {
                        grad[idx_i] = grad[idx_i] + 2 * lambda * bval
                    }
                }
            }
        }

        mom = b1 * mom + (1 - b1) * grad
        vel = b2 * vel + (1 - b2) * (grad:^2)
        step = lr * sqrt(1 - b2^t) / (1 - b1^t)
        delta = delta - step * mom :/ (sqrt(vel) :+ eps)

        if (do_mono) {
            for (i = 1; i <= N; i++) {
                base_i = (i - 1) * p
                for (j = 1; j <= m; j++) {
                    idx_i = base_i + (is_crs ? j : 1 + j)
                    if (delta[idx_i] < 0) delta[idx_i] = 0
                }
            }
        }

        if (verbose & (mod(t, 25) == 0 | t == maxiter)) {
            _convnls_extract(delta, N, m, p, is_crs, alpha, beta)
            mv = 0
            for (k = 1; k <= rows(candidates); k++) {
                i = candidates[k, 1]
                h = candidates[k, 2]
                raw = alpha[i] - alpha[h] + (beta[i,.] - beta[h,.]) * X[h,.]'
                viol = is_concave ? -raw : raw
                if (viol > mv) mv = viol
            }
            printf("{txt}  Adam warm start iter %g: screened max_viol=%9.2e\n",
                   t, mv)
        }
    }

    return(delta)
}


// ============================================================================
// Extract alpha, beta from delta
// ============================================================================
void _convnls_extract(real colvector delta, real scalar N, real scalar m,
                   real scalar p, real scalar is_crs,
                   real colvector alpha, real matrix beta)
{
    real scalar i, base
    alpha = J(N, 1, 0)
    beta  = J(N, m, 0)
    for (i = 1; i <= N; i++) {
        base = (i-1)*p
        if (is_crs) {
            // alpha stays 0; beta occupies positions base+1..base+p
            beta[i,.] = delta[|base+1 \ base+p|]'
        }
        else {
            alpha[i]  = delta[base + 1]
            beta[i,.] = delta[|base+2 \ base+p|]'
        }
    }
}


// ============================================================================
// Find Afriat violations
// Concave: need g_i(x_h) >= g_h(x_h), violation when < -tol
// Convex:  need g_i(x_h) <= g_h(x_h), violation when > +tol
// Returns K x 3: (i, h, |violation|)
// ============================================================================
real matrix _convnls_find_violations(real colvector alpha, real matrix beta,
                                   real matrix X, real scalar N,
                                   real scalar m, real scalar tol,
                                   real scalar is_concave)
{
    real matrix  viol_list
    real colvector v, idx, vals
    real scalar  h, n_viol, max_store, n_take, room

    max_store = min((N*(N-1), 200000))
    viol_list = J(max_store, 3, 0)
    n_viol = 0

    for (h = 1; h <= N; h++) {
        v = alpha :- alpha[h] :+ (beta :- beta[h,.]) * X[h,.]'
        v[h] = 0

        if (is_concave) {
            idx  = selectindex(v :< -tol)
            vals = -v[idx]
        }
        else {
            idx  = selectindex(v :> tol)
            vals = v[idx]
        }

        if (rows(idx) == 0) continue

        room = max_store - n_viol
        if (room <= 0) break
        n_take = min((rows(idx), room))

        viol_list[|n_viol + 1, 1 \ n_viol + n_take, 1|] = idx[1..n_take]
        viol_list[|n_viol + 1, 2 \ n_viol + n_take, 2|] = J(n_take, 1, h)
        viol_list[|n_viol + 1, 3 \ n_viol + n_take, 3|] = vals[1..n_take]
        n_viol = n_viol + n_take

        if (n_viol >= max_store) {
            break
        }
    }

    if (n_viol == 0) return(J(0, 3, 0))
    return(viol_list[1..n_viol, .])
}


// ============================================================================
// Filter violations to only new (i,h) pairs not already in added_pairs
// ============================================================================
real matrix _convnls_filter_new(real matrix violations, real matrix added_pairs)
{
    real scalar k, n_rows, n_keep, i, h
    real colvector keep_idx

    n_rows = rows(violations)
    if (n_rows == 0) return(violations)

    keep_idx = J(n_rows, 1, 0)
    n_keep = 0

    for (k = 1; k <= n_rows; k++) {
        i = violations[k, 1]
        h = violations[k, 2]
        if (added_pairs[i, h] == 0) {
            n_keep++
            keep_idx[n_keep] = k
        }
    }

    if (n_keep == 0) return(J(0, 3, 0))
    return(violations[keep_idx[1..n_keep], .])
}


// ============================================================================
// Mark (i,h) pairs as added
// ============================================================================
void _convnls_mark_added(real matrix violations, real matrix added_pairs)
{
    real scalar k, i, h
    for (k = 1; k <= rows(violations); k++) {
        i = violations[k, 1]
        h = violations[k, 2]
        added_pairs[i, h] = 1
    }
}


// ============================================================================
// Build Afriat constraint rows
// Concave: g_h(x_h) - g_i(x_h) <= 0
// Convex:  g_i(x_h) - g_h(x_h) <= 0
// ============================================================================
real matrix _convnls_build_afriat_rows(real matrix violations, real matrix X,
                                     real scalar N, real scalar m,
                                     real scalar p, real scalar q,
                                     real scalar is_concave,
                                     real scalar is_crs)
{
    real scalar n_rows, k, i, h, j, base_i, base_h
    real matrix rows

    n_rows = rows(violations)
    rows = J(n_rows, N*p + q, 0)

    for (k = 1; k <= n_rows; k++) {
        i = violations[k, 1]
        h = violations[k, 2]
        base_i = (i-1)*p
        base_h = (h-1)*p

        if (is_crs) {
            // No alpha columns; only beta
            if (is_concave) {
                for (j = 1; j <= m; j++) {
                    rows[k, base_h + j] =  X[h, j]
                    rows[k, base_i + j] = -X[h, j]
                }
            }
            else {
                for (j = 1; j <= m; j++) {
                    rows[k, base_i + j] =  X[h, j]
                    rows[k, base_h + j] = -X[h, j]
                }
            }
        }
        else {
            // VRS: include alpha column
            if (is_concave) {
                rows[k, base_h + 1] =  1
                rows[k, base_i + 1] = -1
                for (j = 1; j <= m; j++) {
                    rows[k, base_h + 1 + j] =  X[h, j]
                    rows[k, base_i + 1 + j] = -X[h, j]
                }
            }
            else {
                rows[k, base_i + 1] =  1
                rows[k, base_h + 1] = -1
                for (j = 1; j <= m; j++) {
                    rows[k, base_i + 1 + j] =  X[h, j]
                    rows[k, base_h + 1 + j] = -X[h, j]
                }
            }
        }
    }
    return(rows)
}


// ============================================================================
// Prune inactive Afriat constraints
// Also updates added_pairs to unmark pruned pairs
// ============================================================================
void _convnls_prune(real matrix G_afriat, real matrix G_pairs,
                  real matrix added_pairs,
                  real colvector delta, real scalar prune_tol,
                  real scalar N, real scalar m, real scalar p,
                  real scalar n_pruned)
{
    real colvector slack, keep_idx
    real scalar    n_con, i, n_keep
    real scalar    pi, ph

    n_con = rows(G_afriat)
    if (n_con == 0) {
        n_pruned = 0
        return
    }

    slack = G_afriat * delta

    keep_idx = J(n_con, 1, 0)
    n_keep = 0
    for (i = 1; i <= n_con; i++) {
        if (slack[i] > -prune_tol) {
            n_keep++
            keep_idx[n_keep] = i
        }
        else {
            // Unmark this pair in added_pairs using G_pairs direct lookup
            pi = G_pairs[i, 1]
            ph = G_pairs[i, 2]
            if (pi > 0 & ph > 0) {
                added_pairs[pi, ph] = 0
            }
        }
    }

    n_pruned = n_con - n_keep

    if (n_keep == 0) {
        G_afriat = J(0, cols(G_afriat), 0)
        G_pairs  = J(0, 2, 0)
    }
    else if (n_keep < n_con) {
        G_afriat = G_afriat[keep_idx[1..n_keep], .]
        G_pairs  = G_pairs[keep_idx[1..n_keep], .]
    }
}


// ============================================================================
// ADMM QP Solver with Over-Relaxation
//
// min (1/2) d'Qd - c'd  s.t.  Gd <= 0
// Over-relaxation alpha=1.6, adaptive rho with floor
// ============================================================================
real colvector _convnls_admm_solve(real matrix ZtZ, real colvector Zty,
                                 real matrix G, real colvector delta0,
                                 real scalar rho, real scalar tol,
                                 real scalar maxiter, real scalar verbose)
{
    real scalar    n_vars, n_con, k, i, reg
    real scalar    r_pri, r_dual, eps_pri, eps_dual
    real scalar    alpha_relax, new_rho, prev_r_pri
    real scalar    diverge_count
    real matrix    H, GtG, Gt
    real colvector delta, s, s_old, u, Gdelta, Gdelta_hat, rhs

    n_vars = cols(G)
    n_con  = rows(G)

    if (n_con == 0) {
        return(cholsolve(ZtZ + 1e-8*I(n_vars), Zty))
    }

    alpha_relax = 1.6

    Gt  = G'
    GtG = Gt * G
    reg = 1e-8
    H   = ZtZ + rho * GtG + reg * I(n_vars)

    delta  = delta0
    Gdelta = G * delta
    s      = J(n_con, 1, 0)
    for (i = 1; i <= n_con; i++) {
        s[i] = min((Gdelta[i], 0))
    }
    u = J(n_con, 1, 0)

    r_pri  = .
    r_dual = .
    prev_r_pri = .
    diverge_count = 0

    for (k = 1; k <= maxiter; k++) {

        rhs   = Zty + rho * Gt * (s - u)
        delta = cholsolve(H, rhs)

        if (hasmissing(delta)) {
            delta = lusolve(H, rhs)
            if (hasmissing(delta)) {
                printf("{err}ADMM: solve failed at iter %g\n", k)
                return(delta0)
            }
        }

        Gdelta = G * delta
        Gdelta_hat = alpha_relax * Gdelta + (1 - alpha_relax) * s

        s_old = s
        s     = Gdelta_hat + u
        for (i = 1; i <= n_con; i++) {
            s[i] = min((s[i], 0))
        }

        u = u + Gdelta_hat - s

        // Convergence check every 25 iters
        if (mod(k, 25) == 0 | k <= 5) {
            r_pri  = sqrt(sum((Gdelta - s):^2))
            r_dual = rho * sqrt(sum((Gt * (s - s_old)):^2))

            eps_pri  = sqrt(n_con)  * tol +
                       tol * max((sqrt(sum(Gdelta:^2)), sqrt(sum(s:^2))))
            eps_dual = sqrt(n_vars) * tol +
                       tol * sqrt(sum((rho * Gt * u):^2))

            if (verbose & mod(k, 1000) == 0) {
                printf("{txt}    ADMM %5.0f: r_p=%9.2e r_d=%9.2e rho=%5.2f\n",
                       k, r_pri, r_dual, rho)
            }

            if (r_pri < eps_pri & r_dual < eps_dual) {
                if (verbose) {
                    printf("{txt}    ADMM converged at iter %g (rho=%6.2f)\n", k, rho)
                }
                break
            }

            // Divergence detection: r_pri growing for 3 checks
            if (prev_r_pri < . & r_pri > 10 * prev_r_pri) {
                diverge_count = diverge_count + 1
            }
            else {
                diverge_count = 0
            }
            prev_r_pri = r_pri

            if (diverge_count >= 3) {
                // Reset: reduce over-relaxation and restart u
                if (verbose) {
                    printf("{txt}    ADMM diverging, reset at iter %g\n", k)
                }
                alpha_relax = 1.0  // disable over-relaxation
                u = J(n_con, 1, 0) // reset dual
                s = Gdelta
                for (i = 1; i <= n_con; i++) {
                    s[i] = min((s[i], 0))
                }
                diverge_count = 0
                prev_r_pri = .
            }
        }

        // Adaptive rho every 100 iters; scale u only when rho changes.
        if (mod(k, 100) == 0 & r_pri < . & r_dual < .) {
            if (r_pri > 5 * r_dual & r_dual > 0) {
                new_rho = rho * 2
                if (new_rho != rho) {
                    u   = u :* (rho / new_rho)
                    rho = new_rho
                    H   = ZtZ + rho * GtG + reg * I(n_vars)
                }
            }
            else if (r_dual > 5 * r_pri & r_pri > 0) {
                new_rho = max((rho / 2, 0.1))
                if (new_rho != rho) {
                    u   = u :* (rho / new_rho)
                    rho = new_rho
                    H   = ZtZ + rho * GtG + reg * I(n_vars)
                }
            }
        }
    }

    if (k > maxiter) {
        if (verbose) {
            printf("{txt}    ADMM maxiter (%g)\n", maxiter)
        }
    }

    return(delta)
}


// ============================================================================
// Sort matrix descending by column
// ============================================================================
real matrix _convnls_sort_desc(real matrix M, real scalar col_idx)
{
    real colvector idx
    idx = order(M[., col_idx], -1)
    return(M[idx, .])
}


// ============================================================================
// vcstrat=group: keep only the worst violation per observation i
// O(K) implementation using per-i tracking vectors
// ============================================================================
real matrix _convnls_filter_group(real matrix violations, real scalar N)
{
    real colvector best_viol, best_k
    real scalar    k, i
    real colvector keep_rows

    if (rows(violations) == 0) return(violations)

    best_viol = J(N, 1, 0)
    best_k    = J(N, 1, 0)

    for (k = 1; k <= rows(violations); k++) {
        i = violations[k, 1]
        if (violations[k, 3] > best_viol[i]) {
            best_viol[i] = violations[k, 3]
            best_k[i]    = k
        }
    }

    keep_rows = selectindex(best_k :> 0)
    if (rows(keep_rows) == 0) return(J(0, 3, 0))
    return(violations[best_k[keep_rows], .])
}


// ============================================================================
// Percentile helper: p-th percentile of a colvector (p in [0,100])
// ============================================================================
real scalar _convnls_pctile(real colvector v, real scalar p)
{
    real colvector vs
    real scalar    n, idx

    vs = sort(select(v, v :< .), 1)
    n  = rows(vs)
    if (n == 0) return(0)
    idx = ceil(n * p / 100)
    if (idx < 1) idx = 1
    if (idx > n) idx = n
    return(vs[idx])
}


// ============================================================================
// init=ele: Elementary Afriat initialization
// Sort observations by first scaled input; add 2(n-1) adjacent-pair constraints
// Returns (i, h, 0) triplets; marks pairs in added_pairs
// ============================================================================
real matrix _convnls_init_ele(real matrix X, real scalar N,
                            real matrix added_pairs)
{
    real colvector sort_idx
    real matrix    pairs
    real scalar    k, ia, ib

    sort_idx = order(X[., 1], 1)
    pairs    = J(2*(N-1), 3, 0)

    for (k = 1; k <= N-1; k++) {
        ia = sort_idx[k]
        ib = sort_idx[k+1]
        pairs[2*k-1, .] = (ib, ia, 0)
        pairs[2*k,   .] = (ia, ib, 0)
        added_pairs[ib, ia] = 1
        added_pairs[ia, ib] = 1
    }
    return(pairs)
}


// ============================================================================
// init=sweet: Sweet Spot initialization
// Extends ele with pairs within sweetpct-th percentile of Euclidean distance
// in scaled X space; marks pairs in added_pairs
// ============================================================================
real matrix _convnls_init_sweet(real matrix X, real scalar N, real scalar sweetpct,
                              real matrix added_pairs)
{
    real matrix    pairs, extra, D
    real colvector di
    real scalar    i, j, pd, n_extra, max_extra

    // Start with ele pairs
    pairs = _convnls_init_ele(X, N, added_pairs)

    // Pairwise Euclidean distances in scaled X space
    D = J(N, N, .)
    for (i = 1; i <= N; i++) {
        for (j = i+1; j <= N; j++) {
            D[i,j] = sqrt(sum((X[i,.] - X[j,.]):^2))
            D[j,i] = D[i,j]
        }
    }

    // Pre-allocate buffer (worst case N*(N-1) extra pairs)
    max_extra = N * (N-1)
    extra     = J(max_extra, 3, 0)
    n_extra   = 0

    for (i = 1; i <= N; i++) {
        di    = D[i, .]'
        di[i] = .
        pd    = _convnls_pctile(di, sweetpct)
        if (pd == 0) continue

        for (j = 1; j <= N; j++) {
            if (j == i) continue
            if (D[i,j] > pd) continue
            // Add (i,j) if not already tracked
            if (added_pairs[i,j] == 0) {
                n_extra++
                extra[n_extra,.] = (i, j, 0)
                added_pairs[i,j] = 1
            }
            // Also add reverse (j,i)
            if (added_pairs[j,i] == 0) {
                n_extra++
                extra[n_extra,.] = (j, i, 0)
                added_pairs[j,i] = 1
            }
        }
    }

    if (n_extra > 0) return(pairs \ extra[1..n_extra, .])
    return(pairs)
}

end
