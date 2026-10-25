*! version 0.4.2-rawdiag2  06oct2026
*! stoned - pyStoNED 0.7.3 compatibility with raw invalid-result diagnostics
*!
*! Kuosmanen & Kortelainen (2012) two-stage StoNED estimator.
*! Stage 2: decomposes CNLS residuals into noise (v) and inefficiency (u)
*! default public MOM/QLE and ADDI/MULT conventions; implementation(legacy) uses MM/PSL,
*! then computes firm-level efficiency scores via the
*! Battese & Coelli (1988) conditional expectation E[exp(-u)|z],
*! where z reconstructs the composed error from centered CNLS residuals.
*!
*! SYNTAX: stoned resid_var [if] [in] [, options]
*!
*! Author: Choonjoo Lee (sarang90@korea.kr)

program define _stoned_fit, eclass sortpreserve
    version 16.0
    syntax varname(max=1 numeric) [if] [in], ///
        [FRONTier(string)]   ///
        [METHod(string)]     ///
        [GENerate(name)]     ///
        [INEFFiciency(name)]  ///
        [IMPLEMENTation(string)] [CET(string)] [OUTcome(varname numeric)] ///
        [PYTHONExec(string)] ///
        [INVALID(string)] ///
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

    if "`implementation'" == "" local implementation "pystoned"
    local implementation = lower("`implementation'")
    if !inlist("`implementation'", "pystoned", "legacy") {
        di as error "implementation() must be pystoned or legacy"
        exit 198
    }
    if "`invalid'"=="" local invalid refuse
    local invalid=lower("`invalid'")
    if !inlist("`invalid'","refuse","raw") {
        di as error "invalid() must be refuse or raw"
        exit 198
    }
    if "`implementation'"=="legacy" & "`invalid'"!="refuse" {
        di as error "invalid(raw) is available only with implementation(pystoned)"
        exit 198
    }
    if "`method'" == "" local method "mm"
    local method = lower("`method'")
    if "`method'" == "mom" local method "mm"
    if ("`implementation'" == "pystoned" & !inlist("`method'", "mm", "qle")) | ///
       ("`implementation'" == "legacy" & !inlist("`method'", "mm", "psl")) {
        di as error "pystoned supports method(mm|qle); legacy supports method(mm|psl)"
        exit 198
    }
    if "`implementation'" == "pystoned" {
        if "`cet'" == "" local cet "addi"
        local cet = lower("`cet'")
        if !inlist("`cet'", "addi", "mult") exit 198
        if "`cet'" == "addi" {
            if "`outcome'" == "" {
                di as error "pyStoNED ADDI scoring requires outcome(original_y); use cet(mult) for MULT."
                exit 198
            }
            quietly count if `touse' & (missing(`outcome') | `outcome' == 0)
            if r(N) {
                di as error "outcome() must be nonmissing and nonzero on the entire residual sample."
                exit 498
            }
        }
        else if "`outcome'" != "" {
            di as error "outcome() is only used with cet(addi)"
            exit 198
        }
    }
    else if "`cet'`outcome'" != "" {
        di as error "legacy returns BC scores; cet()/outcome() are not legacy options"
        exit 198
    }

    if "`pythonexec'" != "" & "`method'" != "qle" {
        di as error "pythonexec() is only used by method(qle)"
        exit 198
    }
    if "`generate'" == "" local generate "_stoned_te"
    local outputs `generate' `inefficiency'
    local unique_outputs : list uniq outputs
    local n_outputs : word count `outputs'
    local n_unique : word count `unique_outputs'
    if `n_outputs' != `n_unique' {
        di as error "generate() and inefficiency() must have distinct names"
        exit 198
    }
    local protected `resid_var' `touse' `outcome'
    local collisions : list outputs & protected
    if "`collisions'" != "" {
        di as error "output names conflict with the residual variable: `collisions'"
        exit 198
    }
    foreach v of local outputs {
        capture confirm variable `v', exact
        if _rc confirm new variable `v', exact
    }

    quietly count if `touse'
    if r(N) == 0 error 2000
    local N = r(N)
    if `N' < 3 {
        di as error "At least three nonmissing residuals are required."
        exit 2001
    }

    // Validate provenance before destructive output writes. Characteristics survive
    // save/use, rename and full clonevar; an unmarked residual is never certified.
    local source_status "external_unverified"
    local stage1_converged = .
    local stage1_N = .
    local stage1_frontier ""
    local stage1_rts ""
    local stage1_version ""
    local has_source = 0
    foreach field in schema version converged inner_status frontier rts n stop_reason signature {
        local source_`field' : char `resid_var'[cnls_`field']
        if `"`source_`field''"' != "" local has_source = 1
    }
    if `has_source' {
        if `"`source_schema'"' != "1" | `"`source_version'"' == "" | ///
            !inlist(`"`source_converged'"', "0", "1", "2") | ///
            !inlist(`"`source_frontier'"', "production", "cost") | ///
            !inlist(`"`source_rts'"', "vrs", "crs") | ///
            `"`source_signature'"' == "" | `"`source_stop_reason'"' == "" | ///
            missing(real(`"`source_n'"')) | real(`"`source_n'"') <= 0 {
            di as error "Incomplete or unsupported CNLS residual provenance; refit the first stage."
            exit 498
        }
        if !inlist(`"`source_converged'"', "1", "2") | `"`source_inner_status'"' != "1" {
            di as error "CNLS residuals are not certified (converged=`source_converged')."
            di as error "Refit CNLS to convergence before running stoned."
            exit 498
        }
        if `"`source_frontier'"' != "`frontier'" {
            di as error "frontier() conflicts with the CNLS residual's recorded frontier."
            exit 498
        }
        tempvar signature_resid
        quietly gen double `signature_resid' = `resid_var'
        quietly _datasignature `signature_resid', nonames
        local current_signature "`r(datasignature)'"
        drop `signature_resid'
        quietly count if !missing(`resid_var')
        if `"`current_signature'"' != `"`source_signature'"' | r(N) != real(`"`source_n'"') {
            di as error "CNLS residual values or observation count have changed; refit the first stage."
            exit 498
        }
        local stage1_converged = real(`"`source_converged'"')
        local stage1_N = real(`"`source_n'"')
        local stage1_frontier "`source_frontier'"
        local stage1_rts "`source_rts'"
        local stage1_version "`source_version'"
        local source_status = cond(`stage1_converged'==1, "cnls_strict", "cnls_practical")
    }
    else {
        // Compatibility guard for older CNLS versions which leave no persistent tag.
        if "`e(cmd)'" == "convnls" & "`e(resid_var)'" == "`resid_var'" {
            if !inlist(e(converged),1,2) {
                di as error "The preceding CNLS fit did not certify these residuals."
                exit 498
            }
            if "`e(frontier)'" != "`frontier'" {
                di as error "frontier() conflicts with the preceding CNLS fit."
                exit 498
            }
            local source_status "legacy_unverified"
            local stage1_converged = e(converged)
            local stage1_N = e(N)
            local stage1_frontier "`e(frontier)'"
            local stage1_rts "`e(rts)'"
        }
        di as text "Note: residual provenance is unverified; check first-stage convergence and frontier."
    }

    // Roll back dataset changes if any later calculation fails or is interrupted.
    preserve
    foreach v of local outputs {
        capture confirm variable `v', exact
        if !_rc drop `v'
    }

    // ================================================================
    // 1. Moments of CNLS residuals (needed for MM; used as PSL start)
    // ================================================================
    quietly summ `resid_var' if `touse', detail
    local eps_mean  = r(mean)
    local m2        = r(Var)
    local skew      = r(skewness)
    local sigma_eps = r(sd)
    local m3        = `skew' * (`sigma_eps')^3
    if "`implementation'" == "pystoned" {
        tempvar centered2 centered3
        quietly gen double `centered2' = (`resid_var'-`eps_mean')^2 if `touse'
        quietly gen double `centered3' = (`resid_var'-`eps_mean')^3 if `touse'
        quietly summarize `centered2' if `touse', meanonly
        local m2 = r(mean)
        quietly summarize `centered3' if `touse', meanonly
        local m3 = r(mean)
    }
    local m3_raw = `m3'
    local m3_used = `m3'
    if missing(`m2', `m3') | `m2' <= 0 {
        di as error "Residual moments must be finite with strictly positive variance."
        exit 498
    }

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
    local mm_valid = 1

    if (sign(`m3') != `exp_sign' & `m3' != 0) {
        local wrong_skew = 1
        di as text _n "{bf:Note:} Residual skewness has the wrong sign " ///
            "for `frontier' frontier."
        di as text "  Boundary handling follows implementation(`implementation')."
    }

    if `wrong_skew' {
        local sigma_u_mm = 0
        local sigma_v_mm = sqrt(`m2')
        local m3_used = 0
        if "`implementation'" == "pystoned" & `is_cost' {
            local m3_used = .00001
            local sigma_u_mm = (`m3_used'/`c_mom')^(1/3)
            local sigma_v_mm = sqrt(`m2'-`sigma_u_mm'^2*(1-2/_pi))
            local mm_valid = !missing(`sigma_v_mm') & `sigma_v_mm'>0
        }
    }
    else {
        local sigma_u_mm = (`exp_sign' * `m3' / `c_mom')^(1/3)
        local sigma_v2   = `m2' - `sigma_u_mm'^2 * (1 - 2/_pi)
        local sigma_v_mm = sqrt(`sigma_v2')
        if missing(`sigma_u_mm', `sigma_v_mm') | `sigma_v2' <= 0 {
            local mm_valid = 0
            // Zero noise is not the zero-inefficiency boundary.
            // PSL may still find an interior solution; invalid MM cannot be a fallback.
        }
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
    local psl_boundary_selected = 0
    local psl_ll_gain = .
    local psl_boundary_tol = 1e-8
    local psl_selection "mm_fallback"

    if "`method'" == "psl" {
        // Call-local result scalars must not overwrite users' named scalars.
        tempname su_psl sv_psl conv_psl boundary_psl gain_psl zero_noise_psl
        quietly mata: _stoned_psl("`resid_var'", "`touse'", ///
                                   `is_cost',               ///
                                   `sigma_u_mm', `sigma_v_mm', ///
                                   "`su_psl'", "`sv_psl'", "`conv_psl'", ///
                                   "`boundary_psl'", "`gain_psl'", "`zero_noise_psl'", `psl_boundary_tol')
        if scalar(`zero_noise_psl') {
            di as error "A better zero-noise likelihood limit was found; this boundary is not supported."
            di as error "No PSL efficiency scores were saved."
            exit 498
        }
        local sigma_u_psl   = scalar(`su_psl')
        local sigma_v_psl   = scalar(`sv_psl')
        local psl_converged = scalar(`conv_psl')
        local psl_boundary_selected = scalar(`boundary_psl')
        local psl_ll_gain = scalar(`gain_psl')
        if `psl_converged' {
            local psl_selection = cond(`psl_boundary_selected', "normal_boundary", "interior")
        }
        if !`psl_converged' & `mm_valid' {
            di as text _n "{bf:Note:} PSL optimiser did not converge. " ///
                "Falling back to MOM estimates."
        }
        if `psl_boundary_selected' {
            di as text _n "{bf:Note:} PSL selected sigma_u = 0 after comparing the normal-boundary likelihood."
            di as text "  No converged interior candidate improved mean log-likelihood by more than 1e-8."
        }
        else if `psl_converged' & `sigma_u_psl'/`sigma_v_psl' < 1e-4 {
            di as text _n "{bf:Note:} PSL converged near sigma_u = 0 " ///
                "(positive interior candidate)."
            di as text "  The candidate improved on the normal boundary; inspect e(psl_ll_gain)."
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

    local qle_converged = .
    local qle_lambda = .
    if "`method'" == "qle" {
        tempname qsu qsv qconv qlam
        quietly findfile stoned_py073_qle.py
        local qle_bridge `"`r(fn)'"'
        if "`pythonexec'" == "" local pythonexec "python3"
        // Shell arguments are single-quoted. Refuse embedded quotes/newlines.
        foreach arg in pythonexec qle_bridge {
            if strpos(`"``arg''"', "'") | strpos(`"``arg''"', char(10)) {
                di as error "QLE executable and bridge paths must not contain quotes or newlines."
                exit 198
            }
        }
        tempfile qle_input qle_output
        tempvar bridge_eps
        quietly gen double `bridge_eps' = `resid_var' if `touse'
        format `bridge_eps' %24.17g
        quietly export delimited `bridge_eps' if `touse' using "`qle_input'", replace datafmt
        quietly shell '`pythonexec'' '`qle_bridge'' '`qle_input'' '`qle_output'' `is_cost'
        capture confirm file "`qle_output'"
        if _rc {
            di as error "QLE requires pythonexec(path) with pyStoNED 0.7.3/NumPy 2.0.2/SciPy 1.13.1."
            di as error "Reference/runtime errors are not replaced by native PSL or MM."
            exit 498
        }
        tempname qle_file
        file open `qle_file' using "`qle_output'", read text
        file read `qle_file' qle_header
        if "`qle_header'"!="STONED_QLE_RAW_V1" {
            file close `qle_file'
            di as error "Incomplete or incompatible public QLE response."
            exit 498
        }
        foreach value in qsu qsv qconv qlam {
            file read `qle_file' qle_line
            local `value'_wire "`qle_line'"
            scalar ``value'' = real("`qle_line'")
        }
        file read `qle_file' qle_line
        local qle_status_code "`qle_line'"
        file read `qle_file' qle_line
        local qle_nfev "`qle_line'"
        file read `qle_file' qle_line
        local qle_message "`qle_line'"
        file read `qle_file' qle_line
        local qsu_repr "`qle_line'"
        file read `qle_file' qle_line
        local qsv_repr "`qle_line'"
        file read `qle_file' qle_line
        local qlam_repr "`qle_line'"
        file read `qle_file' qle_footer
        file read `qle_file' qle_extra
        local qle_eof=r(eof)
        file close `qle_file'
        local qle_protocol_ok = ("`qle_footer'"=="END_STONED_QLE_RAW_V1" & `qle_eof')
        if !inlist(scalar(`qconv'),0,1) local qle_protocol_ok=0
        if missing(real("`qle_status_code'"),real("`qle_nfev'")) local qle_protocol_ok=0
        if real("`qle_status_code'")!=floor(real("`qle_status_code'")) local qle_protocol_ok=0
        if real("`qle_nfev'")<1 | real("`qle_nfev'")!=floor(real("`qle_nfev'")) local qle_protocol_ok=0
        foreach value in qsu qsv qlam {
            if missing(scalar(``value'')) {
                if !inlist("``value'_wire'","nan","inf","-inf") local qle_protocol_ok=0
                if "``value'_wire'"!="``value'_repr'" local qle_protocol_ok=0
            }
            else if scalar(``value'')!=real("``value'_repr'") local qle_protocol_ok=0
        }
        if !`qle_protocol_ok' {
            di as error "Malformed public QLE response; raw mode cannot replace a failed bridge."
            exit 498
        }
        if (missing(scalar(`qsu'), scalar(`qsv'), scalar(`qconv'), scalar(`qlam'))) & "`invalid'"!="raw" exit 498
        local sigma_u = scalar(`qsu')
        local sigma_v = scalar(`qsv')
        local qle_converged = scalar(`qconv')
        local qle_lambda = scalar(`qlam')
        if !`qle_converged' & "`invalid'"!="raw" {
            di as error "QLE optimizer did not converge; no fallback is used."
            exit 498
        }
    }

    local effective_method "`method'"
    if "`method'" == "psl" & !`psl_converged' local effective_method "mm"
    local parameter_valid = !missing(`sigma_u',`sigma_v') & `sigma_u'>=0 & `sigma_v'>0
    local invalid_reason "none"
    if missing(`sigma_u') | missing(`sigma_v') local invalid_reason "nonfinite_public_parameter"
    else if `sigma_u'<0 local invalid_reason "negative_sigma_u"
    else if `sigma_v'<=0 local invalid_reason "nonpositive_sigma_v"
    else if "`effective_method'"=="mm" & !`mm_valid' local invalid_reason "invalid_moments_decomposition"
    else if "`method'"=="qle" & !`qle_converged' local invalid_reason "optimizer_not_converged"
    local raw_sigma_u_text : display %24.17g `sigma_u'
    local raw_sigma_v_text : display %24.17g `sigma_v'
    if "`method'"=="mm" & missing(`sigma_v') local raw_sigma_v_text "nan"
    if ("`effective_method'" == "mm" & !`mm_valid') | ///
        !`parameter_valid' | ("`method'"=="qle" & !`qle_converged') {
        if "`implementation'"=="pystoned" & "`invalid'"=="raw" {
            di as text "Reference returned an invalid result; invalid(raw) preserves diagnostics and does not certify its score."
        }
        else {
        di as error "No supported variance decomposition: positive finite noise variance is required."
        di as error "The zero-noise boundary is not implemented; no efficiency scores were saved."
        exit 498
        }
    }
    local score_status "interior"
    if `sigma_u' == 0 local score_status "normal_boundary"

    // ================================================================
    // 4. Display header
    // ================================================================
    local method_label = cond("`effective_method'"=="mm", ///
        "Method of Moments", cond("`method'"=="qle", "pyStoNED QLE", "Pseudo-Likelihood"))
    di as text ""
    di as text "{hline 65}"
    di as text "StoNED: Stage 2 — `method_label' " ///
        "(Kuosmanen & Kortelainen 2012)"
    di as text "{hline 65}"
    di as text "Residual variable: " as result "`resid_var'"
    di as text "Method:            " as result upper("`effective_method'")
    di as text "Frontier type:     " as result ///
        cond(`is_cost', "Cost (z = epsilon + E(u))", ///
                        "Production (z = epsilon - E(u))")
    di as text "Observations:      " as result %6.0f `N'
    di as text "{hline 65}"
    di as text "Results:"
    di as text "{hline 40}"
    di as text "sigma_u (`effective_method'):  " as result %10.6f `sigma_u'
    di as text "sigma_v (`effective_method'):  " as result %10.6f `sigma_v'
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
    //    continued fraction to avoid underflow and subtractive cancellation.
    // ================================================================
    quietly gen double `generate' = . if `touse'
    if "`inefficiency'" != "" {
        quietly gen double `inefficiency' = . if `touse'
    }

    // SCORE_WRITE: stable marker for test-only post-write rollback injection.
    tempname score_valid
    scalar `score_valid' = 1
    if "`implementation'" == "legacy" {
        mata: _stoned_store_scores("`resid_var'", "`touse'", "`generate'", ///
            "`inefficiency'", `sigma_u', `sigma_v', `exp_sign')
    }
    else {
        local rawmode = ("`invalid'"=="raw")
        mata: _stoned_py073_scores("`resid_var'", "`touse'", "`outcome'", ///
            "`generate'", "`inefficiency'", "`cet'", `sigma_u', `sigma_v', ///
            `exp_sign', "`score_valid'", `rawmode')
        if !scalar(`score_valid') {
            quietly count if `touse' & missing(`generate')
            if r(N)>0 {
                local score_status "reference_nonfinite_score"
                di as text "Reference calculation is nonfinite; raw score cells are Stata missing."
            }
            else {
                local score_status "reference_negative_u"
                di as text "Warning: pyStoNED 0.7.3 formula produced negative conditional u; see e(score_valid)."
            }
        }
    }

    label variable `generate' "StoNED efficiency score (`effective_method')"
    if "`inefficiency'" != "" {
        label variable `inefficiency' "StoNED E(u|z) (`effective_method')"
    }
    if "`invalid'"=="raw" {
        label variable `generate' "pyStoNED raw reference score; check e(estimate_valid)"
        if "`inefficiency'"!="" label variable `inefficiency' "pyStoNED raw conditional quantity; check e(estimate_valid)"
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
    // Store the stage-two sample without inventing a coefficient vector.
    ereturn post, esample(`touse') properties("")
    ereturn local implementation "`implementation'"
    ereturn local reference_version = cond("`implementation'"=="pystoned", "0.7.3", "RC3-dev5")
    ereturn local score_definition = cond("`implementation'"=="legacy", "bc", "pystoned_`cet'")
    ereturn local cet "`cet'"
    ereturn local outcome "`outcome'"
    ereturn scalar score_valid = scalar(`score_valid')
    ereturn scalar parameter_valid = `parameter_valid'
    ereturn scalar estimate_valid = `parameter_valid' & scalar(`score_valid') & ///
        ("`method'"!="qle" | `qle_converged'==1)
    ereturn local invalid_policy "`invalid'"
    ereturn local invalid_reason "`invalid_reason'"
    if "`invalid'"=="raw" {
        ereturn local raw_sigma_u_text "`raw_sigma_u_text'"
        ereturn local raw_sigma_v_text "`raw_sigma_v_text'"
    }
    if "`method'"=="qle" {
        ereturn local raw_sigma_u_text "`qsu_repr'"
        ereturn local raw_sigma_v_text "`qsv_repr'"
        ereturn local raw_lambda_text "`qlam_repr'"
        ereturn scalar qle_optimizer_success = `qle_converged'
        ereturn scalar qle_optimizer_status = real("`qle_status_code'")
        ereturn scalar qle_optimizer_nfev = real("`qle_nfev'")
        ereturn local qle_optimizer_message "`qle_message'"
    }
    ereturn scalar m3_raw = `m3_raw'
    ereturn scalar m3_used = `m3_used'
    if "`method'" == "qle" ereturn scalar qle_converged = `qle_converged'
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
        ereturn scalar psl_boundary_selected = `psl_boundary_selected'
        ereturn scalar psl_ll_gain = `psl_ll_gain'
        ereturn scalar psl_boundary_tol = `psl_boundary_tol'
        ereturn local psl_selection "`psl_selection'"
    }
    ereturn local source_status "`source_status'"
    ereturn scalar stage1_converged = `stage1_converged'
    ereturn scalar stage1_N = `stage1_N'
    ereturn local stage1_frontier "`stage1_frontier'"
    ereturn local stage1_rts "`stage1_rts'"
    ereturn local stage1_version "`stage1_version'"
    ereturn scalar mm_valid = `mm_valid'
    ereturn local effective_method "`effective_method'"
    ereturn local score_status "`score_status'"
    ereturn local  method    "`method'"
    ereturn local  resid_var "`resid_var'"
    ereturn local  frontier  "`frontier'"
    ereturn local  te_var    "`generate'"
    if "`inefficiency'" != "" {
        ereturn local u_var   "`inefficiency'"
    }
    ereturn local  cmd       "stoned"
    restore, not
end


// ================================================================
// Mata: PSL negative log-likelihood and optimiser
// ================================================================
mata:

// For t >= 8, q(t) = phi(t)/Phi(-t) - t, evaluated without subtraction.
real scalar _stoned_tail_excess(real scalar t)
{
    real scalar k, tail
    tail = 0
    for (k=100; k>=2; k--) tail = k/(t+tail)
    return(1/(t+tail))
}

// Conditional N(mu, ss^2) truncated below at zero: (E(exp(-u)), E(u)).
// Negative tails use Mills ratios; the mixed-tail identity avoids ss^2 overflow.
real rowvector _stoned_conditional(real scalar mu, real scalar ss)
{
    real scalar a, b, t, q, logte, u
    if (hasmissing((mu,ss)) | ss<=0) return((.,.))
    a = mu/ss
    if (missing(a)) return((.,.))
    b = a-ss
    if (a < -8) {
        t = -a
        q = _stoned_tail_excess(t)
        u = ss*q
        logte = -ln(1+(ss+_stoned_tail_excess(t+ss)-q)/(t+q))
    }
    else {
        u = mu + ss*exp(lnnormalden(a)-lnnormal(a))
        if (b < -8) {
            t = -b
            logte = -.5*a^2 -.5*ln(2*pi()) - ///
                ln(t+_stoned_tail_excess(t)) - lnnormal(a)
        }
        else logte = -mu + ss^2/2 + lnnormal(b)-lnnormal(a)
    }
    if (hasmissing((logte,u)) | logte>1e-12 | u<0) return((.,.))
    return((exp(min((logte,0))),u))
}

void _stoned_store_scores(string scalar resid, string scalar touse,
    string scalar te, string scalar uvar, real scalar su, real scalar sv,
    real scalar sgn)
{
    real colvector eps
    real matrix scores
    real scalar scale, ratio, ss, mu, i
    eps = st_data(.,resid,touse)
    scores = J(rows(eps),2,0)
    if (su==0) scores[,1] = J(rows(eps),1,1)
    else {
        scale = max((su,sv))
        ratio = (su/scale)/sqrt((su/scale)^2+(sv/scale)^2)
        ss = sv*ratio
        for (i=1;i<=rows(eps);i++) {
            mu = (sgn*eps[i]+su*sqrt(2/pi()))*ratio^2
            scores[i,] = _stoned_conditional(mu,ss)
        }
    }
    if (hasmissing(scores)) {
        errprintf("Conditional scores could not be evaluated as finite numbers.\n")
        _error(498)
    }
    st_store(.,te,touse,scores[,1])
    if (uvar!="") st_store(.,uvar,touse,scores[,2])
}

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
    if (hasmissing((su,sv)) | su<=0 | sv<=0) {
        fval = .
        return
    }
    sig = max((su,sv))*sqrt((su/max((su,sv)))^2+(sv/max((su,sv)))^2)
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
    lnPhi = lnnormal(arg)
    // Mata sum() ignores missing entries; reject instead of dropping observations.
    if (hasmissing(lnphi) | hasmissing(lnPhi) | hasmissing((sig,lam))) {
        fval = .
        return
    }

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
                 real   scalar sv_mm,
                 string scalar su_name,
                 string scalar sv_name,
                 string scalar conv_name,
                 string scalar boundary_name,
                 string scalar gain_name,
                 string scalar zero_noise_name,
                 real scalar boundary_tol)
{
    external real matrix  _psl_eps
    external real scalar  _psl_is_cost

    real rowvector theta_hat
    real scalar    rc, su_ret, sv_ret, conv_ret
    real scalar eps_scale, rmse, boundary_ret, gain_ret, f0
    real scalar zero_ret, c, cm, su_zero, f_zero, f_selected
    real colvector t
    transmorphic   S

    // set globals for evaluator
    _psl_eps     = st_data(., varname, tousevar)
    _psl_is_cost = is_cost

    // Optimize in RMSE units: positive rescaling then gives the same problem.
    // Do not center eps: the zero-u candidate has sv^2 = mean(eps^2).
    eps_scale = max(abs(_psl_eps))
    rmse = eps_scale*sqrt(mean((_psl_eps:/eps_scale):^2))
    _psl_eps = _psl_eps:/rmse
    boundary_ret = 0
    gain_ret = .
    zero_ret = 0
    f0 = rows(_psl_eps)*(.5*log(2*pi()) + .5)

    real scalar sv_ols
    sv_ols = sqrt(variance(vec(_psl_eps)))   // start scale, not the boundary MLE

    // Two starting points: MOM estimates (primary) and neutral/small-lambda (fallback)
    real rowvector th_mm, th_neut
    if (!hasmissing((su_mm,sv_mm)) & su_mm/rmse > 1e-12 & sv_mm/rmse > 1e-12) {
        th_mm = (log(su_mm/rmse), log(sv_mm/rmse))
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
    optimize_init_conv_ptol(S, 1e-9)
    optimize_init_conv_vtol(S, 1e-9)
    optimize_init_conv_nrtol(S, 1e-10)
    optimize_init_conv_maxiter(S, 500)
    optimize_init_verbose(S, 0)
    rc = _optimize(S)

    // ---- attempt 2: neutral simplex start for a flat near-zero-u likelihood ----
    // Always run; keep whichever converged solution has lower neg-LL (= higher LL)
    real scalar rc2
    transmorphic S2
    S2 = optimize_init()
    optimize_init_evaluator(S2, &_stoned_psl_negll())
    optimize_init_evaluatortype(S2, "d0")
    optimize_init_params(S2, th_neut)
    optimize_init_which(S2, "min")
    optimize_init_technique(S2, "nm")
    optimize_init_nmsimplexdeltas(S2, (.1,.1))
    optimize_init_conv_ptol(S2, 1e-9)
    optimize_init_conv_vtol(S2, 1e-9)
    optimize_init_conv_nrtol(S2, 1e-10)
    optimize_init_conv_maxiter(S2, 500)
    optimize_init_verbose(S2, 0)
    rc2 = _optimize(S2)

    // A lower-noise start checks the distinct mode missed by the first two starts.
    real scalar rc3
    transmorphic S3
    S3 = optimize_init()
    optimize_init_evaluator(S3, &_stoned_psl_negll())
    optimize_init_evaluatortype(S3, "d0")
    optimize_init_params(S3, (log(1.5),log(.1)))
    optimize_init_which(S3, "min")
    optimize_init_technique(S3, "bfgs")
    optimize_init_conv_ptol(S3, 1e-9)
    optimize_init_conv_vtol(S3, 1e-9)
    optimize_init_conv_nrtol(S3, 1e-10)
    optimize_init_conv_maxiter(S3, 500)
    optimize_init_verbose(S3, 0)
    rc3 = _optimize(S3)

    // Select: prefer converged solution with lower neg-LL (= higher log-likelihood)
    real scalar ok1, ok2, ok3
    ok1 = (rc == 0 & optimize_result_converged(S)==1 & optimize_result_value(S)<.)
    ok2 = (rc2 == 0 & optimize_result_converged(S2)==1 & optimize_result_value(S2)<.)
    ok3 = (rc3 == 0 & optimize_result_converged(S3)==1 & optimize_result_value(S3)<.)
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
    if (ok3) {
        if (!(ok1 | ok2)) S = S3
        else if (optimize_result_value(S3) < optimize_result_value(S)) S = S3
    }

    if (ok1 | ok2 | ok3) {
        theta_hat = optimize_result_params(S)
        su_ret    = rmse*exp(theta_hat[1])
        sv_ret    = rmse*exp(theta_hat[2])
        conv_ret  = 1
        // Compare only after an interior optimization actually converged.
        // A finite boundary alone must not mask failure of all starts.
        gain_ret = (f0 - optimize_result_value(S))/rows(_psl_eps)
        f_selected = optimize_result_value(S)
        if (gain_ret <= boundary_tol) {
            su_ret = 0
            sv_ret = rmse
            boundary_ret = 1
            f_selected = f0
        }
        // Reject when the unsupported zero-noise limit beats our selected candidate.
        // At a support endpoint use the limit approached from inside the support.
        if (_psl_is_cost) t = _psl_eps
        else t = -_psl_eps
        c = sqrt(2/pi())
        cm = c*mean(t)
        su_zero = max(((cm+sqrt(cm^2+4))/2, -min(t)/c))
        f_zero = rows(t)*(log(su_zero)-log(2)+.5*log(2*pi())) + .5*sum((t:/su_zero:+c):^2)
        if ((f_selected-f_zero)/rows(t) > boundary_tol) zero_ret = 1
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
    st_numscalar(su_name,   su_ret)
    st_numscalar(sv_name,   sv_ret)
    st_numscalar(conv_name, conv_ret)
    st_numscalar(boundary_name, boundary_ret)
    st_numscalar(gain_name, gain_ret)
    st_numscalar(zero_noise_name, zero_ret)
}

end


mata:
// Literal published 0.7.3 score convention, including its 1e-6 denominator.
// This is compatibility, not a stable-tail correction or BC score.
void _stoned_py073_scores(string scalar rv, string scalar tv, string scalar yv,
    string scalar svname, string scalar uv, string scalar cet,
    real scalar su, real scalar sv, real scalar sign, string scalar validname,
    real scalar rawmode)
{
    real colvector e, a, u, score, y
    real scalar sig, ss
    e=st_data(.,rv,tv):+sign*su*sqrt(2/pi())
    sig=sqrt(su^2+sv^2)
    ss=su*sv/sig
    a=e:*su/(sv*sig)
    u=ss:*(normalden(a):/(1:-normal(-sign:*a):+1e-6):+sign:*a)
    if (cet=="addi") {
        y=st_data(.,yv,tv)
        score=(y:+sign:*u):/y
    }
    else score=exp(sign:*u)
    if (hasmissing((u,score))) {
        if (!rawmode) _error(498)
        st_numscalar(validname,0)
    }
    st_store(.,svname,tv,score)
    if (uv!="") st_store(.,uv,tv,u)
    if (!hasmissing((u,score))) st_numscalar(validname,min(u)>=0)
}

end
