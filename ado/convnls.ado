*! version 0.5.3-dev11  29sep2026
*! convnls - Convex Nonparametric Least Squares (StoNED Stage 1)
*!
*! Kuosmanen (2008, 2012) CNLS with Cutting Plane + ADMM.
*!
*! SYNTAX: convnls depvar [if] [in], INputs(varlist) [options]
*!
*! Author: Choonjoo Lee (sarang90@korea.kr)

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
        [RHOPOLicy(string)] ///
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

    if `n_obs' < 2 {
        di as error "at least two observations are required"
        exit 2001
    }
    if missing(`tolerance') | `tolerance' <= 0 {
        di as error "tolerance() must be finite and positive"
        exit 198
    }
    if missing(`rho') | (`rho' != -1 & `rho' <= 0) {
        di as error "rho() must be -1 (automatic) or finite and positive"
        exit 198
    }
    if missing(`maxiter') | missing(`maxcut') | `maxiter' < 1 | `maxcut' < 1 {
        di as error "maxiter() and maxcut() must be positive integers"
        exit 198
    }
    if missing(`prune') | `prune' < 0 {
        di as error "prune() must be finite and nonnegative"
        exit 198
    }
    if !inrange(`sweetpct', 0, 100) {
        di as error "sweetpct() must be between 0 and 100"
        exit 198
    }

    if "`rhopolicy'" == "" local rhopolicy "normalized"
    local rhopolicy = lower("`rhopolicy'")
    if !inlist("`rhopolicy'", "normalized", "raw") {
        di as error "rhopolicy() must be normalized or raw"
        exit 198
    }

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

    if "`vce'" == "" local vce "none"
    local vce = lower("`vce'")
    if !inlist("`vce'", "none", "ols", "robust") {
        di as error "vce() must be {bf:none}, {bf:ols}, or {bf:robust}"
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
    if missing(`adamiter') | `adamiter' < 0 {
        di as error "adamiter() must be nonnegative"
        exit 198
    }
    if missing(`adamlr') | `adamlr' <= 0 {
        di as error "adamlr() must be positive"
        exit 198
    }
    if missing(`adamlambda') | `adamlambda' <= 0 {
        di as error "adamlambda() must be positive"
        exit 198
    }
    if missing(`adamscreen') | `adamscreen' < 0 {
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

    if missing(`addper') | `addper' < 0 {
        di as error "addper() must be nonnegative"
        exit 198
    }

    if "`generate'" == "" local generate "_convnls_resid"
    if "`fitted'" == ""   local fitted "_convnls_fitted"

    // Validate the entire write set before touching existing variables.
    local outputs `generate' `fitted' _convnls_alpha
    forvalues j = 1/`m' {
        local outputs `outputs' _convnls_beta`j'
    }
    local unique_outputs : list uniq outputs
    local n_outputs : word count `outputs'
    local n_unique : word count `unique_outputs'
    if `n_outputs' != `n_unique' {
        di as error "output names must be distinct from each other and the coefficient names"
        exit 198
    }
    local protected `depvar' `inputs' `context' `touse'
    local collisions : list outputs & protected
    if "`collisions'" != "" {
        di as error "output names conflict with input variables: `collisions'"
        exit 198
    }
    foreach v of local outputs {
        capture confirm variable `v', exact
        if _rc confirm new variable `v', exact
    }

    foreach v of varlist `depvar' `inputs' `context' {
        quietly count if missing(`v') & `touse'
        if r(N) > 0 {
            di as error "Variable `v' has `r(N)' missing values"
            exit 459
        }
    }

    // A design check is necessary, but does not certify statistical identification
    // or sampling inference for the jointly estimated nonparametric frontier.
    if `q' > 0 {
        local response_context : list depvar & context
        if "`response_context'" != "" {
            di as error "the dependent variable cannot be included in context()"
            exit 198
        }
        mata: _convnls_context_check("`inputs'", "`context'", "`touse'")
    }

    // Stata restores this dataset automatically on error or Break. Commit only
    // at normal program completion. Exact matching avoids abbreviation drops.
    preserve
    foreach v of local outputs {
        capture confirm variable `v', exact
        if !_rc drop `v'
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
        di as text "Covariance diagnostic: " as result upper("`vce'")
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
        "`verbose'", "`rhopolicy'")

    // Save Mata r() results before any r-class command overwrites them.
    tempname mata_cp_history mata_admm_info
    matrix `mata_admm_info' = r(admm_info)
    local mata_max_mono = r(max_mono)
    local mata_inner_iters = r(total_inner_iters)
    local mata_inner_solves = r(inner_solves)
    local mata_warm_starts = r(warm_starts)
    local mata_final_rho = r(final_rho)
    local mata_constraint_bytes = r(constraint_bytes)
    local mata_dense_bytes = r(dense_constraint_bytes)
    tempname mata_factor_stats
    matrix `mata_factor_stats' = r(factor_stats)
    local mata_stop_reason "`r(stop_reason)'"
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

    // Undefined and exactly zero skewness carry no directional sign evidence.
    if missing(`skew') {
        di as text ""
        di as text "Residual skewness is undefined; no sign comparison is available."
    }
    else if `skew' == 0 {
        di as text ""
        di as text "Residual skewness = 0; no directional skewness is indicated."
    }
    else if `is_concave' {
        if `skew' > 0 {
            di as text ""
            di as text "{bf:Note:} Residual skewness > 0 (" ///
                as result %9.3g `skew' as text ")."
            di as text "  Production frontier expects negative skewness."
        }
        else {
            di as text ""
            di as text "Residual skewness < 0 (" as result %9.3g `skew' ///
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
            di as text "Residual skewness > 0 (" as result %9.3g `skew' ///
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

    // Persist stage-one status on the residual itself, independently of e().
    // Double copy normalizes storage type; _datasignature is order independent.
    tempvar signature_resid
    quietly gen double `signature_resid' = `generate'
    quietly _datasignature `signature_resid', nonames
    local residual_signature "`r(datasignature)'"
    drop `signature_resid'
    char `generate'[cnls_schema] "1"
    char `generate'[cnls_version] "0.5.3-dev11"
    char `generate'[cnls_converged] "`mata_converged'"
    char `generate'[cnls_inner_status] "`= `mata_admm_info'[1,1]'"
    char `generate'[cnls_frontier] "`frontier'"
    char `generate'[cnls_rts] "`rts'"
    char `generate'[cnls_n] "`n_obs'"
    char `generate'[cnls_stop_reason] "`mata_stop_reason'"
    char `generate'[cnls_signature] "`residual_signature'"

    ereturn clear
    if `q' > 0 {
        // Point estimates only: do not expose an unvalidated covariance as e(V).
        ereturn post `context_b', esample(`touse') properties(b)
    }
    else {
        // No parametric coefficient vector: post only the actual estimation sample.
        ereturn post, esample(`touse') properties("")
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
    ereturn scalar max_mono = `mata_max_mono'
    ereturn scalar inner_status = `mata_admm_info'[1,1]
    ereturn scalar inner_iter = `mata_admm_info'[1,2]
    ereturn scalar inner_primal = `mata_admm_info'[1,3]
    ereturn scalar inner_dual = `mata_admm_info'[1,4]
    ereturn scalar inner_eps_primal = `mata_admm_info'[1,5]
    ereturn scalar inner_eps_dual = `mata_admm_info'[1,6]
    ereturn scalar inner_stationarity = `mata_admm_info'[1,7]
    ereturn scalar inner_eps_stationarity = `mata_admm_info'[1,8]
    ereturn scalar inner_complementarity = `mata_admm_info'[1,9]
    ereturn scalar inner_eps_complementarity = `mata_admm_info'[1,10]
    ereturn scalar inner_tolerance = `mata_admm_info'[1,11]
    ereturn scalar total_inner_iters = `mata_inner_iters'
    ereturn scalar inner_solves = `mata_inner_solves'
    ereturn scalar warm_starts = `mata_warm_starts'
    ereturn local rho_policy "`rhopolicy'"
    ereturn scalar rho_adapt_limit = cond("`rhopolicy'" == "raw", 0, 5000)
    ereturn scalar final_rho = `mata_final_rho'
    ereturn scalar constraint_bytes = `mata_constraint_bytes'
    ereturn scalar dense_constraint_bytes = `mata_dense_bytes'
    ereturn local engine "structured"
    ereturn scalar factorizations = `mata_factor_stats'[1,1]
    ereturn scalar chol_solves = `mata_factor_stats'[1,2]
    ereturn scalar lu_fallbacks = `mata_factor_stats'[1,3]
    ereturn scalar rho_updates = `mata_factor_stats'[1,4]
    ereturn local stop_reason "`mata_stop_reason'"
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
    ereturn local vce          "none"
    ereturn local residual_schema "1"
    ereturn local residual_signature "`residual_signature'"
    ereturn local version "0.5.3-dev11"
    ereturn local cmd          "convnls"
    if `q' > 0 {
        ereturn local inference_status "not_validated"
        ereturn local context_diagnostic "`vce'"
        ereturn local context_diagnostic_status "not_requested"
        if "`vce'" != "none" {
            ereturn local context_diagnostic_status "not_certified"
            if inlist(`mata_converged',1,2) & `mata_admm_info'[1,1]==1 {
                ereturn matrix V_context_diagnostic = `context_V'
                ereturn matrix V_robust_diagnostic = `context_V_robust'
                ereturn local context_diagnostic_status "legacy_formula_only"
            }
        }
        di as text ""
        di as text "Context-variable point estimates (joint additive CNLS):"
        matrix list e(b), noheader
        di as text "SEs, tests, and confidence intervals are not validated; e(V) is not posted."
        if "`vce'" != "none" {
            di as text "Requested legacy covariance diagnostic: `vce'; status: `e(context_diagnostic_status)'"
        }
    }
    restore, not
end


// ============================================================================
// MATA ENGINE
// ============================================================================
version 16.0
mata:
mata set matastrict on

// Relative SVD threshold after unit normalization. This detects numerical
// linear dependence; passing it is not a proof of nonlinear identification.
real scalar _convnls_context_rank(real matrix A)
{
    real colvector s
    s = svdsv(A)
    if (max(s)==0) return(0)
    return(sum(s :> 1e-10*max(s)))
}

void _convnls_context_check(string scalar inputs, string scalar context,
                          string scalar touse)
{
    real matrix X, C
    real scalar j, scale, rx
    X = st_data(., inputs, touse)
    C = st_data(., context, touse)
    for (j=1; j<=cols(X); j++) {
        X[.,j] = X[.,j] :- mean(X[.,j])
        scale = sqrt(sum(X[.,j]:^2))
        if (scale>0) X[.,j] = X[.,j]/scale
    }
    for (j=1; j<=cols(C); j++) {
        C[.,j] = C[.,j] :- mean(C[.,j])
        scale = sqrt(sum(C[.,j]:^2))
        if (scale>0) C[.,j] = C[.,j]/scale
    }
    if (_convnls_context_rank(C)<cols(C)) {
        errprintf("context() contains constant or numerically dependent columns\n")
        _error(459)
    }
    rx = _convnls_context_rank(X)
    if (_convnls_context_rank((X,C))<rx+cols(C)) {
        errprintf("context() is numerically aliased with an affine function of inputs()\n")
        _error(459)
    }
}

// Compact two-block Afriat operator. Context columns are identically zero.
struct _convnls_operator {
    real scalar dense, n, p, q, sign, nvars, ncons, raw_rho
    real matrix G, Gt, Z, pairs, zh, zhi, pi, ph
    real colvector oi, oh, ui, uh, mono
}

// Continuation belongs to one fixed QP and its returned primal iterate.
// The caller must invalidate it whenever Q, c, or the ordered constraints change.
struct _convnls_admm_state {
    real scalar valid, reused, nvars, ncons, rho, iteration, pending_rho, raw_rho
    real scalar alpha, prev_primal, diverge_count, primal, dual, eps_primal, eps_dual
    real colvector s, u
}

// Balance distance from the existing stopping thresholds, not raw residual units.
// These ratios guide the penalty only; all four certification checks stay intact.
real scalar _convnls_next_rho(real scalar rho, real scalar primal, real scalar dual,
                            real scalar eps_primal, real scalar eps_dual)
{
    real scalar rp, rd
    if (hasmissing((rho,primal,dual,eps_primal,eps_dual))) return(rho)
    if (rho<=0 | primal<0 | dual<0 | eps_primal<=0 | eps_dual<=0) return(rho)
    rp=primal/eps_primal
    rd=dual/eps_dual
    if (hasmissing((rp,rd))) return(rho)
    // Bounds must not reverse the requested direction for an explicit rho()
    // that starts outside them. Zero residuals are valid, including both zero.
    if (rp/5 > rd & rho < 1e6) return(min((rho*2,1e6)))
    if (rd/5 > rp & rho > 1e-6) return(max((rho/2,1e-6)))
    return(rho)
}

// Raw is the exact dev4 residual-balancing policy, including its 0.1 floor.
real scalar _convnls_penalty(real scalar rho, real scalar primal, real scalar dual,
                           real scalar eps_primal, real scalar eps_dual,
                           real scalar raw_rho)
{
    if (!raw_rho) return(_convnls_next_rho(rho,primal,dual,eps_primal,eps_dual))
    if (primal > 5*dual & dual > 0) return(rho*2)
    if (dual > 5*primal & primal > 0) return(max((rho/2,0.1)))
    return(rho)
}

struct _convnls_operator scalar _convnls_pair_operator(real matrix X,
    real matrix pairs, real scalar q, real scalar is_crs,
    real scalar is_concave, real scalar do_mono)
{
    struct _convnls_operator scalar op
    real scalar i, j, k, m, a
    real colvector ids
    op.dense=0
    op.raw_rho=0
    op.n=rows(X)
    m=cols(X)
    op.p=m+(!is_crs)
    op.q=q
    op.sign=is_concave ? -1 : 1
    op.nvars=op.n*op.p+q
    op.Z=is_crs ? X : (J(op.n,1,1),X)
    op.pairs=pairs
    a=rows(pairs)
    op.mono=J(do_mono*op.n*m,1,0)
    k=0
    if (do_mono) for (i=1;i<=op.n;i++) for (j=1;j<=m;j++) {
        k++
        op.mono[k]=(i-1)*op.p+(!is_crs)+j
    }
    op.ncons=a+rows(op.mono)
    op.G=op.Gt=J(0,0,0)
    op.zh=op.zhi=J(0,op.p,0)
    op.pi=op.ph=J(0,2,0)
    op.oi=op.oh=op.ui=op.uh=J(0,1,0)
    if (a) {
        op.zh=op.Z[pairs[,2],.]
        op.oi=order(pairs,1)
        op.oh=order(pairs,2)
        ids=pairs[op.oi,1]
        op.pi=panelsetup(ids,1)
        op.ui=ids[op.pi[,1]]
        ids=pairs[op.oh,2]
        op.ph=panelsetup(ids,1)
        op.uh=ids[op.ph[,1]]
        op.zhi=op.zh[op.oi,.]
    }
    return(op)
}

real colvector _convnls_op_times(struct _convnls_operator scalar op,
                                real colvector d)
{
    real matrix D
    real colvector out
    real scalar a
    if (op.dense) return(op.G*d)
    if (hasmissing(d)) return(J(op.ncons,1,.))
    a=rows(op.pairs)
    out=J(op.ncons,1,0)
    if (a) {
        D=colshape(d[|1\op.n*op.p|],op.p)
        out[|1\a|]=op.sign*rowsum(op.zh:*(D[op.pairs[,1],.]-D[op.pairs[,2],.]),1)
    }
    if (rows(op.mono)) out[|a+1\op.ncons|]=-d[op.mono]
    return(out)
}

real colvector _convnls_op_transpose(struct _convnls_operator scalar op,
                                    real colvector v)
{
    real matrix D
    real colvector out
    real scalar a
    if (op.dense) return(op.Gt*v)
    if (hasmissing(v)) return(J(op.nvars,1,.))
    a=rows(op.pairs)
    D=J(op.n,op.p,0)
    if (a) {
        D[op.ui,.]=panelsum(op.zhi,v[op.oi],op.pi)
        D[op.uh,.]=D[op.uh,.]-panelsum(v[op.oh],op.ph):*op.Z[op.uh,.]
    }
    out=(op.sign*vec(D') \ J(op.q,1,0))
    if (rows(op.mono)) out[op.mono]=out[op.mono]-v[|a+1\op.ncons|]
    return(out)
}

real matrix _convnls_op_gram(struct _convnls_operator scalar op)
{
    real matrix out, zz
    real rowvector ii, hh
    real scalar k, a
    if (op.dense) return(op.Gt*op.G)
    out=J(op.nvars,op.nvars,0)
    a=rows(op.pairs)
    for (k=1;k<=a;k++) {
        ii=(op.pairs[k,1]-1)*op.p:+(1..op.p)
        hh=(op.pairs[k,2]-1)*op.p:+(1..op.p)
        zz=op.zh[k,.]'*op.zh[k,.]
        out[ii,ii]=out[ii,ii]+zz
        out[hh,hh]=out[hh,hh]+zz
        out[ii,hh]=out[ii,hh]-zz
        out[hh,ii]=out[hh,ii]-zz
    }
    for (k=1;k<=rows(op.mono);k++) out[op.mono[k],op.mono[k]]=out[op.mono[k],op.mono[k]]+1
    return(out)
}

real scalar _convnls_op_bytes(struct _convnls_operator scalar op)
{
    // Matrix payload only, excluding allocator overhead and temporary workspaces.
    return(8*(length(op.Z)+length(op.pairs)+length(op.zh)+length(op.zhi)+
              length(op.pi)+length(op.ph)+length(op.oi)+length(op.oh)+
              length(op.ui)+length(op.uh)+length(op.mono)))
}

void _convnls_prune_pairs(real matrix pairs, real matrix added,
    real matrix X, real colvector d, real scalar q, real scalar is_crs,
    real scalar is_concave, real scalar tolerance, real scalar removed)
{
    struct _convnls_operator scalar op
    real colvector slack, keep
    real scalar k
    op=_convnls_pair_operator(X,pairs,q,is_crs,is_concave,0)
    slack=_convnls_op_times(op,d)
    keep=selectindex(slack:> -tolerance)
    for (k=1;k<=rows(pairs);k++) {
        if (slack[k]<= -tolerance) added[pairs[k,1],pairs[k,2]]=0
    }
    removed=rows(pairs)-rows(keep)
    pairs=rows(keep) ? pairs[keep,.] : J(0,2,0)
}



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
                    string scalar verbose_str, string scalar rho_policy)
{
    real matrix    X, X_raw, C, C_raw, C_center, ZtZ
    real matrix    beta_curr, violations, cp_history
    real colvector y, y_raw, obs_idx, Zty, delta, alpha_curr
    real colvector resid, resid_scaled, fitted, x_scale, c_mean, c_scale
    real colvector theta_scaled, theta_orig
    real matrix    context_V, context_V_robust, XtX_inv, meat
    real rowvector zi, ci
    real scalar    N, m, q, p, p_local, verbose, i, j, base, cbase
    real scalar    cp_iter, max_viol, n_active_afriat
    real scalar    resid_idx, fitted_idx, alpha_idx, beta_j_idx
    real scalar    n_mono_viol, rho, inner_tol
    real scalar    y_mean, y_sd, n_pruned, n_new_unique
    real scalar near_converge_count, converged, n_viol, max_mono
    real scalar total_inner_iters, inner_solves, warm_starts
    struct _convnls_admm_state scalar admm_state
    real rowvector admm_info
    real rowvector factor_stats
    string scalar stop_reason
    // Constraint tracking: N x N matrix, 1 if (i,h) already added
    real matrix    added_pairs
    real scalar    is_crs
    real matrix    G_pairs
    real scalar peak_operator_bytes, peak_dense_bytes
    struct _convnls_operator scalar op
    real matrix    init_pairs

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
    delta = _convnls_init_ols(y, X, C, N, m, q, p_local, is_crs, do_mono)
    if (verbose) printf("{txt}OLS init done, rho=%6.2f\n", rho)

    // ---- Initialize constraint tracking ----
    added_pairs = J(N, N, 0)

    // ---- Cutting plane: solve first, then certify the returned iterate ----
    G_pairs = J(0, 2, 0)
    n_active_afriat = 0
    near_converge_count = 0
    converged = 0
    stop_reason = "maxcut"
    total_inner_iters = 0
    inner_solves = 0
    warm_starts = 0
    admm_state.valid = 0
    admm_info = J(1, 15, .)
    factor_stats = J(1, 4, 0)
    peak_operator_bytes=peak_dense_bytes=0
    max_viol = .
    n_new_unique = 0

    if (init_str != "none") {
        if (init_str == "ele") {
            init_pairs = _convnls_init_ele(X, N, added_pairs)
        }
        else if (init_str == "sweet") {
            init_pairs = _convnls_init_sweet(X, N, sweetpct_val, added_pairs)
        }
        else {
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
            G_pairs = init_pairs[., 1..2]
        }
    }

    printf("{txt}\nCutting Plane:\n")
    printf("{txt}  %4s  %7s  %7s  %7s  %10s\n",
           "Iter", "Active", "Violat", "NewAdd", "Max Viol")
    cp_history = J(0, 5, .)

    for (cp_iter = 1; cp_iter <= maxcut; cp_iter++) {
        // Any pruning is followed by a solve; an old certificate cannot survive it.
        if (prune_tol > 0 & rows(G_pairs) > 0 & cp_iter > 3) {
            _convnls_prune_pairs(G_pairs, added_pairs, X, delta, q, is_crs,
                                is_concave, prune_tol, n_pruned)
            if (n_pruned > 0) admm_state.valid = 0
        }
        n_active_afriat = rows(G_pairs)
        op=_convnls_pair_operator(X,G_pairs,q,is_crs,is_concave,do_mono)
        op.raw_rho=(rho_policy=="raw")
        peak_operator_bytes=max((peak_operator_bytes,_convnls_op_bytes(op)))
        peak_dense_bytes=max((peak_dense_bytes,8*op.ncons*op.nvars))
        inner_tol = tol * 0.1
        if (max_viol < . & max_viol >= 30*tol & n_new_unique > 0) {
            inner_tol = max((0.01*max_viol, tol))
        }
        delta = _convnls_admm_resume(ZtZ, Zty, op, delta, rho,
                                  inner_tol, maxiter, verbose, admm_info, admm_state)
        total_inner_iters = total_inner_iters + admm_info[2]
        inner_solves++
        warm_starts = warm_starts + admm_state.reused
        factor_stats = factor_stats + admm_info[12..15]
        _convnls_extract(delta, N, m, p_local, is_crs, alpha_curr, beta_curr)
        violations = _convnls_separate(alpha_curr, beta_curr, X, N, tol,
                                      is_concave, added_pairs, vcstrat_str,
                                      addper, max_viol, n_viol)
        max_mono = do_mono ? max((0, -min(beta_curr))) : 0

        // A loose intermediate solve is never a convergence certificate.
        if (max_viol < 30*tol & max_mono <= tol & inner_tol > tol & admm_info[1] >= 0) {
            delta = _convnls_admm_resume(ZtZ, Zty, op, delta, rho,
                                      tol*0.1, maxiter, verbose, admm_info, admm_state)
            total_inner_iters = total_inner_iters + admm_info[2]
            inner_solves++
            warm_starts = warm_starts + admm_state.reused
            factor_stats = factor_stats + admm_info[12..15]
            _convnls_extract(delta, N, m, p_local, is_crs, alpha_curr, beta_curr)
            violations = _convnls_separate(alpha_curr, beta_curr, X, N, tol,
                                          is_concave, added_pairs, vcstrat_str,
                                          addper, max_viol, n_viol)
            max_mono = do_mono ? max((0, -min(beta_curr))) : 0
        }

        if (admm_info[1] == 1 & admm_info[11] <= tol & max_mono <= tol & max_viol < 30*tol) {
            near_converge_count++
            if (max_viol <= tol) {
                converged = 1
                stop_reason = "converged"
            }
            else if (near_converge_count >= 2) {
                converged = 2
                stop_reason = "practical"
            }
        }
        else near_converge_count = 0

        n_new_unique = rows(violations)
        if (converged > 0 | cp_iter == maxcut | admm_info[1] < 0) n_new_unique = 0
        cp_history = cp_history \ (cp_iter, n_active_afriat, n_viol, n_new_unique, max_viol)
        printf("{txt}  %4.0f  %7.0f  %7.0f  %7.0f  %10.2e\n",
               cp_iter, n_active_afriat, n_viol, n_new_unique, max_viol)
        if (admm_info[1] < 0) {
            stop_reason = "inner_failure"
            break
        }
        if (converged > 0 | cp_iter == maxcut) break

        if (n_new_unique > 0) {
            _convnls_mark_added(violations, added_pairs)
            G_pairs = G_pairs \ violations[., 1..2]
            admm_state.valid = 0
        }
        // With no new constraints, retain the complete ADMM state and rho.
    }

    // All diagnostics refer to the last returned iterate, including budget exits.
    if (converged == 0 & admm_info[1] == 0) stop_reason = "inner_maxiter"
    st_numscalar("r(converged)", converged)
    st_numscalar("r(cp_iters)", cp_iter)
    st_numscalar("r(n_active)", n_active_afriat)
    st_numscalar("r(max_viol)", max_viol)
    st_numscalar("r(max_mono)", max_mono)
    st_numscalar("r(total_inner_iters)", total_inner_iters)
    st_numscalar("r(inner_solves)", inner_solves)
    st_numscalar("r(warm_starts)", warm_starts)
    st_numscalar("r(final_rho)", rho)
    st_numscalar("r(constraint_bytes)",peak_operator_bytes)
    st_numscalar("r(dense_constraint_bytes)",peak_dense_bytes)
    st_matrix("r(factor_stats)", factor_stats)
    st_global("r(stop_reason)", stop_reason)
    st_matrix("r(cp_history)", cp_history)
    st_matrix("r(admm_info)", admm_info)
    if (converged == 0) printf("{txt}  Not converged (%s).\n", stop_reason)
    else if (converged == 1) printf("{txt}  Converged: inner optimality and global shape checks passed.\n")
    else printf("{txt}  Practically converged: inner optimality passed; global Afriat violation < 30*tol.\n")

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
                               real scalar is_crs, real scalar do_mono)
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
                delta[base + j] = do_mono ? max((scale_i * xi[j], 0)) : scale_i * xi[j]
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
                delta[base + 1 + j] = do_mono ? max((beta_ols[j], 0)) : beta_ols[j]
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
// Complete separation with bounded candidate storage. Global diagnostics include
// active pairs and violations below tol; candidate selection includes NEW pairs only.
real matrix _convnls_separate(real colvector alpha, real matrix beta,
                             real matrix X, real scalar N, real scalar tol,
                             real scalar is_concave, real matrix added_pairs,
                             string scalar strategy, real scalar cap,
                             real scalar max_viol, real scalar n_viol)
{
    real scalar h, k, i, limit, n_keep, take
    real colvector v, idx, best, best_h, ord
    real matrix candidates, batch
    max_viol = 0
    n_viol = 0
    candidates = J(0, 3, 0)
    best = J(N, 1, 0)
    best_h = J(N, 1, 0)
    limit = cap > 0 ? min((cap, N*(N-1))) : min((200000, N*(N-1)))
    if (strategy == "one") limit = 1
    for (h = 1; h <= N; h++) {
        v = alpha :- alpha[h] :+ (beta :- beta[h,.]) * X[h,.]'
        if (is_concave) v = -v
        v[h] = 0
        max_viol = max((max_viol, max(v)))
        n_viol = n_viol + sum(v :> tol)
        idx = selectindex((v :> tol) :& (added_pairs[.,h] :== 0))
        if (rows(idx) == 0) continue
        if (strategy == "group") {
            for (k = 1; k <= rows(idx); k++) {
                i = idx[k]
                if (v[i] > best[i]) {
                    best[i] = v[i]
                    best_h[i] = h
                }
            }
        }
        else {
            batch = (idx, J(rows(idx), 1, h), v[idx])
            candidates = candidates \ batch
            // Reserve at most two candidate batches before selecting the top K.
            // Unlike the old cap, this NEVER terminates the full scan.
            if (rows(candidates) > 2*limit) {
                ord = order(candidates, (-3, 1, 2))
                candidates = candidates[ord[1..limit], .]
            }
        }
    }
    if (strategy == "group") {
        idx = selectindex(best_h :> 0)
        if (rows(idx) == 0) return(J(0, 3, 0))
        return((idx, best_h[idx], best[idx]))
    }
    n_keep = rows(candidates)
    if (n_keep > limit) {
        ord = order(candidates, (-3, 1, 2))
        take = min((n_keep, limit))
        candidates = candidates[ord[1..take], .]
    }
    return(candidates)
}

// Warm-start callers need only candidates, but still scan the entire sample.
real matrix _convnls_find_violations(real colvector alpha, real matrix beta,
                                   real matrix X, real scalar N,
                                   real scalar m, real scalar tol,
                                   real scalar is_concave)
{
    real matrix empty_pairs
    real scalar max_viol, n_viol
    empty_pairs = J(N, N, 0)
    return(_convnls_separate(alpha, beta, X, N, tol, is_concave,
                            empty_pairs, "all", 0, max_viol, n_viol))
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
// info: status iterations primal dual eps_primal eps_dual stationarity
//       eps_stationarity complementarity eps_complementarity tolerance
//       factorizations chol_solves lu_fallbacks rho_updates
// status: 1 certified, 0 iteration limit, -1 numerical failure.
real colvector _convnls_admm_solve(real matrix ZtZ, real colvector Zty,
                                 real matrix G, real colvector delta0,
                                 real scalar rho, real scalar tol,
                                 real scalar maxiter, real scalar verbose,
                                 real rowvector info)
{
    struct _convnls_operator scalar op
    op.dense=1
    op.raw_rho=0
    op.G=G
    op.Gt=G'
    op.nvars=cols(G)
    op.ncons=rows(G)
    return(_convnls_admm_core(ZtZ,Zty,op,delta0,rho,tol,maxiter,verbose,info))
}

real colvector _convnls_admm_core(real matrix ZtZ, real colvector Zty,
                                 struct _convnls_operator scalar op, real colvector delta0,
                                 real scalar rho, real scalar tol,
                                 real scalar maxiter, real scalar verbose,
                                 real rowvector info)
{
    struct _convnls_admm_state scalar state
    state.valid=0
    return(_convnls_admm_resume(ZtZ,Zty,op,delta0,rho,tol,maxiter,verbose,info,state))
}

// Bound normalized adaptation within each unchanged QP. Continuing the same
// state keeps this clock; rebuilding the working set starts a new clock.
// A constant tail avoids sustained residual oscillation from penalty switching.
// raw deliberately retains the dev4 schedule. This is not an iteration limit.
real scalar _convnls_penalty_due(real scalar iteration, real scalar raw_rho)
{
    return(iteration > 0 & mod(iteration,100)==0 & (raw_rho==1 | iteration<5000))
}

real colvector _convnls_admm_resume(real matrix ZtZ, real colvector Zty,
                                 struct _convnls_operator scalar op, real colvector delta0,
                                 real scalar rho, real scalar tol,
                                 real scalar maxiter, real scalar verbose,
                                 real rowvector info, struct _convnls_admm_state scalar state)
{
    real scalar n_vars, n_con, k, i, reg, r_pri, r_dual, eps_pri, eps_dual
    real scalar alpha_relax, new_rho, prev_r_pri, diverge_count
    real scalar stationarity, eps_stat, comp, eps_comp
    real scalar factor_dirty, factor_ok, offset, iteration, regular_check, raw_rho
    real matrix H, GtG, L, Lt
    real colvector delta, delta_prev, s, s_old, u, Gdelta, Gdelta_hat, rhs, grad, dual_term
    raw_rho=(op.raw_rho==1)
    n_vars = op.nvars
    n_con = op.ncons
    info = (0, 0, ., ., ., ., ., ., ., ., tol, 0, 0, 0, 0)
    state.reused=0
    if (n_con == 0) {
        state.valid=0
        delta = invsym(ZtZ) * Zty
        if (hasmissing(delta)) {
            info[1] = -1
            return(delta0)
        }
        grad = ZtZ * delta
        stationarity = max(abs(grad - Zty))
        eps_stat = tol + tol*max((max(abs(grad)), max(abs(Zty))))
        info[1..11] = ((stationarity <= eps_stat ? 1 : -1), 0, 0, 0, tol, tol,
                stationarity, eps_stat, 0, tol, tol)
        return(delta)
    }

    delta = delta0
    Gdelta = _convnls_op_times(op, delta)
    if (state.valid == 1) {
        if (state.nvars != n_vars | state.ncons != n_con | state.raw_rho != raw_rho) state.valid=0
        else if (rows(state.s)!=n_con | rows(state.u)!=n_con) state.valid=0
        else if (hasmissing(state.s) | hasmissing(state.u)) state.valid=0
    }
    if (state.valid == 1) {
        s=state.s
        u=state.u
        rho=state.rho
        offset=state.iteration
        alpha_relax=state.alpha
        prev_r_pri=state.prev_primal
        diverge_count=state.diverge_count
        state.reused=1
        // Defer a penalty update at a budget boundary until continuation starts.
        if (state.pending_rho & _convnls_penalty_due(state.iteration,raw_rho)) {
            new_rho=_convnls_penalty(rho,state.primal,state.dual,
                                      state.eps_primal,state.eps_dual,raw_rho)
            if (new_rho != rho) {
                u=u:*(rho/new_rho)
                rho=new_rho
                info[15]=info[15]+1
            }
        }
    }
    else {
        s = rowmin((Gdelta, J(n_con, 1, 0)))
        u = J(n_con, 1, 0)
        offset=0
        alpha_relax=1.6
        prev_r_pri=.
        diverge_count=0
    }
    GtG = _convnls_op_gram(op)
    reg = 1e-8
    H = ZtZ + rho * GtG + reg * I(n_vars)
    factor_dirty = 1
    factor_ok = 0
    for (k = 1; k <= maxiter; k++) {
        iteration=offset+k
        delta_prev = delta
        // Proximal stabilization: matched RHS term preserves the unregularized QP.
        rhs = Zty + rho * _convnls_op_transpose(op, s - u) + reg * delta_prev
        // H changes only when rho changes within this restricted solve.
        // A new function call starts a fresh cache for its new working set.
        if (factor_dirty) {
            L = cholesky(H)
            info[12] = info[12] + 1
            factor_ok = !hasmissing(L)
            if (factor_ok) {
                // Match cholsolve()'s default singularity threshold, rather
                // than silently accepting a generalized triangular solution.
                factor_ok = min(diagonal(L)) > 1e-13*mean(abs(diagonal(L)))
            }
            if (factor_ok) Lt = L'
            factor_dirty = 0
        }
        if (factor_ok) {
            delta = solveupper(Lt, solvelower(L, rhs))
            info[13] = info[13] + 1
        }
        else delta = J(n_vars, 1, .)
        if (hasmissing(delta)) {
            info[14] = info[14] + 1
            delta = lusolve(H, rhs)
        }
        if (hasmissing(delta)) {
            state.valid=0
            info[1] = -1
            info[2] = k
            info[3..10] = J(1, 8, .)
            return(delta_prev)
        }
        Gdelta = _convnls_op_times(op, delta)
        Gdelta_hat = alpha_relax * Gdelta + (1-alpha_relax)*s
        s_old = s
        s = rowmin((Gdelta_hat + u, J(n_con, 1, 0)))
        u = u + Gdelta_hat - s

        // Always measure the last iterate, even when maxiter is not a multiple of 25.
        regular_check=(mod(iteration,25)==0 | iteration<=5)
        if (regular_check | k == maxiter) {
            dual_term = rho * _convnls_op_transpose(op, u)
            grad = ZtZ * delta
            r_pri = max(abs(Gdelta - s))
            r_dual = rho * max(abs(_convnls_op_transpose(op, s - s_old)))
            eps_pri = tol + tol * max((max(abs(Gdelta)), max(abs(s))))
            eps_dual = tol + tol * max(abs(dual_term))
            stationarity = max(abs(grad - Zty + dual_term))
            eps_stat = tol + tol * max((max(abs(grad)), max(abs(Zty)), max(abs(dual_term))))
            comp = max(abs(s :* (rho*u)))
            eps_comp = tol + tol * max(abs(s)) * max(abs(rho*u))
            info[1..11] = (0, k, r_pri, r_dual, eps_pri, eps_dual,
                    stationarity, eps_stat, comp, eps_comp, tol)
            if (r_pri <= eps_pri & r_dual <= eps_dual &
                stationarity <= eps_stat & comp <= eps_comp) {
                info[1] = 1
                break
            }
            if (verbose & mod(k,1000) == 0) {
                printf("{txt}    ADMM %g: primal=%9.2e dual=%9.2e stationarity=%9.2e\n",
                       k, r_pri, r_dual, stationarity)
            }
            // Final-only diagnostics must not alter the continuation schedule.
            if (regular_check) {
                if (prev_r_pri < . & r_pri > 10*prev_r_pri) diverge_count++
                else diverge_count = 0
                prev_r_pri = r_pri
            }
            if (diverge_count >= 3) {
                alpha_relax = 1
                u = J(n_con,1,0)
                s = rowmin((Gdelta, J(n_con,1,0)))
                diverge_count = 0
                prev_r_pri = .
            }
        }
        if (_convnls_penalty_due(iteration,raw_rho) & k < maxiter) {
            new_rho = _convnls_penalty(rho,r_pri,r_dual,eps_pri,eps_dual,raw_rho)
            if (new_rho != rho) {
                u = u :* (rho/new_rho)
                rho = new_rho
                H = ZtZ + rho*GtG + reg*I(n_vars)
                factor_dirty = 1
                info[15] = info[15] + 1
            }
        }
    }
    state.valid=1
    state.raw_rho=raw_rho
    state.nvars=n_vars
    state.ncons=n_con
    state.s=s
    state.u=u
    state.rho=rho
    state.iteration=offset+info[2]
    state.alpha=alpha_relax
    state.prev_primal=prev_r_pri
    state.diverge_count=diverge_count
    state.primal=info[3]
    state.dual=info[4]
    state.eps_primal=info[5]
    state.eps_dual=info[6]
    state.pending_rho=(info[1]==0 & _convnls_penalty_due(state.iteration,raw_rho))
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
