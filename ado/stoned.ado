*! version 0.3.2  26jun2026
*! stoned - StoNED Stage 2: Method of Moments or Pseudo-Likelihood
*!
*! Kuosmanen & Kortelainen (2012) two-stage StoNED estimator.
*! Stage 2: decomposes CNLS residuals into noise (v) and inefficiency (u)
*! using method of moments (MM) or pseudo-likelihood (PSL),
*! then computes firm-level efficiency scores via the
*! Battese & Coelli (1988) conditional expectation E[exp(-u)|z],
*! where z reconstructs the composed error from centered CNLS residuals.
*!
*! SYNTAX: stoned resid_var [if] [in] [, options]
*!
*! Author: Choonjoo Lee (bloom.rampike@gmail.com)

program define stoned, eclass sortpreserve
    version 16.0
    syntax varname(max=1 numeric) [if] [in], ///
        [FRONTier(string)]   ///
        [METHod(string)]     ///
        [GENerate(name)]     ///
        [INEFFiciency(name)]  ///
        [VERBose]

    marksample touse
    local resid_var `varlist'

    if "`frontier'" == "" local frontier "production"
    local frontier = lower("`frontier'")
    if !inlist("`frontier'", "production", "cost") {
        di as error "frontier() must be {bf:production} or {bf:cost}"
        exit 198
    }
    local is_cost = ("`frontier'" == "cost")

    if "`method'" == "" local method "mm"
    local method = lower("`method'")
    if !inlist("`method'", "mm", "psl") {
        di as error "method() must be {bf:mm} (method of moments) or {bf:psl} (pseudo-likelihood)"
        exit 198
    }

    if "`generate'" == "" local generate "_stoned_te"
    capture drop `generate'
    if "`inefficiency'" != "" {
        capture drop `inefficiency'
    }

    quietly count if `touse'
    if r(N) == 0 error 2000
    local N = r(N)

    // ================================================================
    // 1. Moments of CNLS residuals (needed for MM; used as PSL start)
    // ================================================================
    quietly summ `resid_var' if `touse', detail
    local eps_mean  = r(mean)
    local m2        = r(Var)
    local skew      = r(skewness)
    local sigma_eps = r(sd)
    local m3        = `skew' * (`sigma_eps')^3

    if "`verbose'" != "" {
        di as text "  Residual mean:     " as result %10.6f `eps_mean'
        di as text "  Residual variance: " as result %10.6f `m2'
        di as text "  Residual skewness: " as result %10.6f `skew'
        di as text "  m3 (third moment): " as result %10.6e `m3'
    }

    // ================================================================
    // 2. Method of Moments (MM)
    //
    //   Half-normal model:  m3 = c_mom * sigma_u^3
    //   c_mom = (4/pi - 1)*sqrt(2/pi) ≈ 0.2180
    //   sigma_u = (|m3| / c_mom)^(1/3)
    //   sigma_v = sqrt(m2 - sigma_u^2*(1 - 2/pi))
    // ================================================================
    local c_mom    = (4/_pi - 1) * sqrt(2/_pi)
    local exp_sign = cond(`is_cost', 1, -1)
    local wrong_skew = 0

    if (sign(`m3') != `exp_sign' & `m3' != 0) {
        local wrong_skew = 1
        di as text _n "{bf:Note:} Residual skewness has the wrong sign " ///
            "for `frontier' frontier."
        di as text "  Cannot identify sigma_u via MOM."
        di as text "  Setting sigma_u = 0 (all variation attributed to noise)."
    }

    if `wrong_skew' {
        local sigma_u_mm = 0
        local sigma_v_mm = sqrt(`m2')
    }
    else {
        local sigma_u_mm = (`exp_sign' * `m3' / `c_mom')^(1/3)
        local sigma_v2   = `m2' - `sigma_u_mm'^2 * (1 - 2/_pi)
        if (`sigma_v2' < 0) {
            di as text _n "{bf:Note:} Implied sigma_v^2 < 0 — setting sigma_v = 0."
            local sigma_v2 = 0
        }
        local sigma_v_mm = sqrt(`sigma_v2')
    }

    // ================================================================
    // 3. PSL: Pseudo-likelihood via Mata optimize()
    //    Maximises log L(su,sv) = n*log(2) - n*log(sigma)
    //                           + sum log_phi(z(su)/sigma)
    //                           + sum lnnormal(±lambda*z(su)/sigma)
    //    Parameterised as (log su, log sv) for unconstrained optimisation.
    //    MOM estimates used as starting values.
    //
    //    CNLS residuals are centered around the estimated conditional mean.
    //    The likelihood is therefore evaluated at the reconstructed composed
    //    error z = eps + E(u) for cost and z = eps - E(u) for production.
    // ================================================================
    local psl_converged = 0
    local sigma_u_psl   = `sigma_u_mm'
    local sigma_v_psl   = `sigma_v_mm'

    if "`method'" == "psl" {
        // Mata writes results to Stata scalars _su_psl, _sv_psl, _conv_psl
        quietly mata: _stoned_psl("`resid_var'", "`touse'", ///
                                   `is_cost',               ///
                                   `sigma_u_mm', `sigma_v_mm')
        local sigma_u_psl   = scalar(_su_psl)
        local sigma_v_psl   = scalar(_sv_psl)
        local psl_converged = scalar(_conv_psl)
        scalar drop _su_psl _sv_psl _conv_psl
        if !`psl_converged' {
            di as text _n "{bf:Note:} PSL optimiser did not converge. " ///
                "Falling back to MOM estimates."
        }
        if `psl_converged' & `sigma_u_psl' < 1e-4 {
            di as text _n "{bf:Note:} PSL converged to sigma_u = 0 " ///
                "(boundary solution)."
            di as text "  CNLS residuals may not retain the composed-error " ///
                "mean structure required by PSL."
            di as text "  Consider using {bf:method(mm)} instead."
        }
    }

    // ---- select active estimates ----
    if "`method'" == "mm" {
        local sigma_u = `sigma_u_mm'
        local sigma_v = `sigma_v_mm'
    }
    else {
        local sigma_u = `sigma_u_psl'
        local sigma_v = `sigma_v_psl'
    }

    // ================================================================
    // 4. Display header
    // ================================================================
    local method_label = cond("`method'"=="mm", ///
        "Method of Moments", "Pseudo-Likelihood")
    di as text ""
    di as text "{hline 65}"
    di as text "StoNED: Stage 2 — `method_label' " ///
        "(Kuosmanen & Kortelainen 2012)"
    di as text "{hline 65}"
    di as text "Residual variable: " as result "`resid_var'"
    di as text "Method:            " as result upper("`method'")
    di as text "Frontier type:     " as result ///
        cond(`is_cost', "Cost (z = epsilon + E(u))", ///
                        "Production (z = epsilon - E(u))")
    di as text "Observations:      " as result %6.0f `N'
    di as text "{hline 65}"
    di as text "Results:"
    di as text "{hline 40}"
    di as text "sigma_u (`method'):  " as result %10.6f `sigma_u'
    di as text "sigma_v (`method'):  " as result %10.6f `sigma_v'
    if (`sigma_v' > 0) {
        di as text "lambda (su/sv):      " as result %10.4f `sigma_u'/`sigma_v'
    }
    local sigma2_tot = `sigma_u'^2 + `sigma_v'^2
    if (`sigma2_tot' > 0) {
        di as text "Share of u (%%):      " as result ///
            %10.1f (`sigma_u'^2 / `sigma2_tot' * 100)
    }
    if "`method'" == "psl" {
        di as text "PSL converged:       " as result `psl_converged'
    }

    // ================================================================
    // 5. Conditional efficiency and inefficiency scores
    //    Battese & Coelli (1988): E[exp(-u_i) | epsilon_i]
    //    Jondrow et al. (1982):   E[u_i | epsilon_i]
    //
    //    sigma^2  = sigma_v^2 + sigma_u^2
    //    sigma*   = sigma_v * sigma_u / sigma
    //    z_i      = epsilon_i + sign * E(u_i), sign=+1 cost, -1 production
    //    mu*_i    = sign * sigma_u^2 * z_i / sigma^2
    //
    //    TE_i = exp(-mu*_i + sigma*^2/2) *
    //           Phi(mu*_i/sigma* - sigma*) / Phi(mu*_i/sigma*)
    //
    //    E[u_i | z_i] = mu*_i + sigma* *
    //           phi(mu*_i/sigma*) / Phi(mu*_i/sigma*)
    //    The inverse Mills ratio is evaluated with a left-tail
    //    approximation to avoid underflow when Phi(.) is tiny.
    // ================================================================
    quietly gen double `generate' = . if `touse'
    if "`inefficiency'" != "" {
        quietly gen double `inefficiency' = . if `touse'
    }

    if `sigma_u' > 0 & `sigma_v' > 0 {
        local sigma2 = `sigma_u'^2 + `sigma_v'^2
        local sigma  = sqrt(`sigma2')
        local sstar  = `sigma_v' * `sigma_u' / `sigma'
        local eu      = `sigma_u' * sqrt(2/_pi)
        tempvar comp_err
        quietly gen double `comp_err' = ///
            `resid_var' + `exp_sign' * `eu' if `touse'

        quietly replace `generate' = ///
            exp( -(`exp_sign' * `sigma_u'^2 * `comp_err' / `sigma2') ///
                 + `sstar'^2 / 2 ) * ///
            normal( (`exp_sign' * `sigma_u'^2 * `comp_err' / `sigma2') ///
                    / `sstar' - `sstar' ) / ///
            normal( (`exp_sign' * `sigma_u'^2 * `comp_err' / `sigma2') ///
                    / `sstar' ) ///
            if `touse'

        if "`inefficiency'" != "" {
            tempvar mu_u arg_u imr_u
            quietly gen double `mu_u' = ///
                (`exp_sign' * `sigma_u'^2 * `comp_err' / `sigma2') ///
                if `touse'
            quietly gen double `arg_u' = `mu_u' / `sstar' if `touse'
            quietly gen double `imr_u' = ///
                normalden(`arg_u') / normal(`arg_u') if `touse' & `arg_u' > -8
            quietly replace `imr_u' = -`arg_u' + 1/(-`arg_u') ///
                if `touse' & `arg_u' <= -8
            quietly replace `inefficiency' = `mu_u' + `sstar' * `imr_u' ///
                if `touse'
            quietly replace `inefficiency' = max(`inefficiency', 0) if `touse'
        }

        quietly replace `generate' = 1 if missing(`generate') & `touse'
        quietly replace `generate' = min(`generate', 1) if `touse'
        quietly replace `generate' = max(`generate', 0) if `touse'
    }
    else {
        quietly replace `generate' = 1 if `touse'
        if "`inefficiency'" != "" {
            quietly replace `inefficiency' = 0 if `touse'
        }
    }

    label variable `generate' "StoNED efficiency score (`method')"
    if "`inefficiency'" != "" {
        label variable `inefficiency' "StoNED E(u|z) (`method')"
    }

    quietly summ `generate' if `touse'
    local te_mean = r(mean)
    local te_sd   = r(sd)
    local te_min  = r(min)
    local te_max  = r(max)

    local u_mean = .
    local u_sd   = .
    local u_min  = .
    local u_max  = .
    if "`inefficiency'" != "" {
        quietly summ `inefficiency' if `touse'
        local u_mean = r(mean)
        local u_sd   = r(sd)
        local u_min  = r(min)
        local u_max  = r(max)
    }

    di as text "{hline 40}"
    di as text "Efficiency scores (`generate'):"
    di as text "  Mean:            " as result %10.6f `te_mean'
    di as text "  Std. Dev.:       " as result %10.6f `te_sd'
    di as text "  Min:             " as result %10.6f `te_min'
    di as text "  Max:             " as result %10.6f `te_max'
    if "`inefficiency'" != "" {
        di as text "{hline 40}"
        di as text "Inefficiency estimates (`inefficiency'):"
        di as text "  Mean:            " as result %10.6f `u_mean'
        di as text "  Std. Dev.:       " as result %10.6f `u_sd'
        di as text "  Min:             " as result %10.6f `u_min'
        di as text "  Max:             " as result %10.6f `u_max'
    }
    di as text "{hline 65}"

    // ================================================================
    // 6. Store e() results
    // ================================================================
    ereturn clear
    ereturn scalar N            = `N'
    ereturn scalar sigma_u      = `sigma_u'
    ereturn scalar sigma_v      = `sigma_v'
    if (`sigma_v' > 0) {
        ereturn scalar lambda   = `sigma_u' / `sigma_v'
    }
    ereturn scalar m2           = `m2'
    ereturn scalar m3           = `m3'
    ereturn scalar te_mean      = `te_mean'
    ereturn scalar te_min       = `te_min'
    ereturn scalar te_max       = `te_max'
    if "`inefficiency'" != "" {
        ereturn scalar u_mean   = `u_mean'
        ereturn scalar u_min    = `u_min'
        ereturn scalar u_max    = `u_max'
    }
    ereturn scalar sigma_u_mm   = `sigma_u_mm'
    ereturn scalar sigma_v_mm   = `sigma_v_mm'
    ereturn scalar eu           = `sigma_u' * sqrt(2/_pi)
    if "`method'" == "psl" {
        ereturn scalar sigma_u_psl  = `sigma_u_psl'
        ereturn scalar sigma_v_psl  = `sigma_v_psl'
        ereturn scalar psl_converged = `psl_converged'
    }
    ereturn local  method    "`method'"
    ereturn local  resid_var "`resid_var'"
    ereturn local  frontier  "`frontier'"
    ereturn local  te_var    "`generate'"
    if "`inefficiency'" != "" {
        ereturn local u_var   "`inefficiency'"
    }
    ereturn local  cmd       "stoned"
end


// ================================================================
// Mata: PSL negative log-likelihood and optimiser
// ================================================================
mata:

// ----------------------------------------------------------------
// Negative log-likelihood evaluator (d0 type for optimize())
//
// optimize() d0 evaluator: sets the scalar objective value fval.
// Data (eps, is_cost) accessed via external Mata globals.
// ----------------------------------------------------------------
void _stoned_psl_negll(transmorphic   S,
                        real rowvector theta,
                        real scalar    fval,
                        real rowvector g,
                        real matrix    H)
{
    external real matrix  _psl_eps
    external real scalar  _psl_is_cost

    real scalar su, sv, sig, lam, n, ll
    real matrix comp, z, lnphi, lnPhi, arg

    su  = exp(theta[1,1])
    sv  = exp(theta[1,2])
    sig = sqrt(su^2 + sv^2)
    lam = su / sv
    n   = rows(_psl_eps)

    if (_psl_is_cost) {
        comp = _psl_eps :+ su * sqrt(2/pi())
    }
    else {
        comp = _psl_eps :- su * sqrt(2/pi())
    }

    z      = comp :/ sig
    // log phi(z): numerically stable form
    lnphi  = -0.5*log(2*pi()) :- 0.5 :* z:^2
    // log Phi: Phi(+lam*z) for cost (eps=v+u, positive skewness)
    //          Phi(-lam*z) for production (eps=v-u, negative skewness)
    if (_psl_is_cost) {
        arg = lam :* z
    }
    else {
        arg = -lam :* z
    }
    lnPhi = log(rowmax((normal(arg), J(n, 1, 1e-300))))

    ll   = n*log(2) - n*log(sig) + sum(lnphi) + sum(lnPhi)
    fval = -ll   // minimize negative log-likelihood
}

// ----------------------------------------------------------------
// Run PSL optimisation; writes results back to Stata scalars
// ----------------------------------------------------------------
void _stoned_psl(string scalar varname,
                 string scalar tousevar,
                 real   scalar is_cost,
                 real   scalar su_mm,
                 real   scalar sv_mm)
{
    external real matrix  _psl_eps
    external real scalar  _psl_is_cost

    real rowvector theta0, theta_hat
    real scalar    rc, su_ret, sv_ret, conv_ret
    pointer        S

    // set globals for evaluator
    _psl_eps     = st_data(., varname, tousevar)
    _psl_is_cost = is_cost

    real scalar n, sv_ols
    n      = rows(_psl_eps)
    sv_ols = sqrt(variance(vec(_psl_eps)))   // sigma_v under H0: u=0

    // Two starting points: MOM estimates (primary) and neutral/small-lambda (fallback)
    real rowvector th_mm, th_neut
    if (su_mm > 1e-12 & sv_mm > 1e-12) {
        th_mm = (log(su_mm), log(sv_mm))
    }
    else {
        th_mm = (log(sv_ols * 0.5), log(sv_ols))
    }
    th_neut = (log(sv_ols * 0.2), log(sv_ols * 0.95))  // small su, sv ≈ RMSE

    // ---- attempt 1: MOM starting values (primary) ----
    S = optimize_init()
    optimize_init_evaluator(S, &_stoned_psl_negll())
    optimize_init_evaluatortype(S, "d0")
    optimize_init_params(S, th_mm)
    optimize_init_which(S, "min")
    optimize_init_technique(S, "bfgs")
    optimize_init_conv_ptol(S, 1e-6)
    optimize_init_conv_vtol(S, 1e-6)
    optimize_init_conv_maxiter(S, 500)
    optimize_init_verbose(S, 0)
    rc = _optimize(S)

    // ---- attempt 2: neutral starting values ----
    // Always run; keep whichever converged solution has lower neg-LL (= higher LL)
    real scalar rc2
    pointer S2
    S2 = optimize_init()
    optimize_init_evaluator(S2, &_stoned_psl_negll())
    optimize_init_evaluatortype(S2, "d0")
    optimize_init_params(S2, th_neut)
    optimize_init_which(S2, "min")
    optimize_init_technique(S2, "bfgs")
    optimize_init_conv_ptol(S2, 1e-6)
    optimize_init_conv_vtol(S2, 1e-6)
    optimize_init_conv_maxiter(S2, 500)
    optimize_init_verbose(S2, 0)
    rc2 = _optimize(S2)

    // Select: prefer converged solution with lower neg-LL (= higher log-likelihood)
    real scalar ok1, ok2
    ok1 = (rc  == 0 | rc  == 1 | rc  == 430)
    ok2 = (rc2 == 0 | rc2 == 1 | rc2 == 430)
    if (ok1 & ok2) {
        // Both converged — pick the one with higher log-likelihood
        if (optimize_result_value(S2) < optimize_result_value(S)) {
            S  = S2
            rc = rc2
        }
    }
    else if (!ok1 & ok2) {
        S  = S2
        rc = rc2
    }
    // else: keep attempt 1 result (or both failed)

    if (rc == 0 | rc == 1 | rc == 430) {
        // rc=0: full convergence; rc=1: max iter; rc=430: not achieved (usable)
        theta_hat = optimize_result_params(S)
        su_ret    = exp(theta_hat[1])
        sv_ret    = exp(theta_hat[2])
        conv_ret  = (rc == 0 ? 1 : 2)   // 1=full, 2=partial
    }
    else {
        su_ret   = su_mm
        sv_ret   = sv_mm
        conv_ret = 0
    }

    // clean up globals
    _psl_eps     = J(0,0,.)
    _psl_is_cost = .

    // write results to Stata scalars for retrieval
    st_numscalar("_su_psl",   su_ret)
    st_numscalar("_sv_psl",   sv_ret)
    st_numscalar("_conv_psl", conv_ret)
}

end
