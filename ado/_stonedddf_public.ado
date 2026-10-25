program define _stonedddf_public, rclass
    version 16
    syntax , GENerate(name) METHod(string) [PYTHONExec(string) ALIGNED REPLACE INVALID(string)]
    if "`invalid'"=="" local invalid refuse
    if !inlist("`invalid'","refuse","raw") exit 198
    tempvar sample eps weighted score rawu
    quietly _stonedddf_state, sample(`sample')
    local protected `r(protected)'
    local d `r(distance)'
    local outputs `r(outputs)'
    local coef `r(coefficient)'
    if "`coef'"=="" {
        di as error "public DDF scoring requires same-fit gamma coefficients; refit with coefficients(prefix)"
        exit 459
    }
    local destinations `generate'_score `generate'_raw_u
    local overlap : list destinations & protected
    if "`overlap'"!="" exit 198
    foreach v of local destinations {
        confirm name `v'
        if strlen("`v'")>32 exit 198
        if "`replace'"=="" confirm new variable `v'
    }
    quietly gen double `weighted'=0 if `sample'
    local j=0
    foreach y of local outputs {
        local ++j
        quietly replace `weighted'=`weighted'+`coef'_gamma`j'*`y' if `sample'
    }
    quietly count if `sample' & (missing(`weighted') | `weighted'<=0)
    if r(N) {
        di as error "same-fit gamma-weighted output must be finite and positive"
        exit 498
    }
    quietly summarize `weighted' if `sample', meanonly
    tempname denominator center su sv mu N conv um sm neg
    scalar `denominator'=r(min)
    quietly summarize `d' if `sample', meanonly
    scalar `center'=r(mean)
    quietly gen double `eps'=`d' if `sample'
    if "`aligned'"!="" quietly replace `eps'=-(`d'-scalar(`center')) if `sample'
    local pyopt
    if "`method'"=="qle" & `"`pythonexec'"'!="" local pyopt pythonexec(`"`pythonexec'"')
    tempname held
    _estimates hold `held', restore copy
    capture noisily _stoned_fit `eps' if `sample', implementation(pystoned) ///
        method(`method') frontier(production) cet(addi) outcome(`weighted') ///
        generate(`score') inefficiency(`rawu') `pyopt' invalid(`invalid')
    local code=_rc
    if !`code' {
        scalar `su'=e(sigma_u)
        scalar `sv'=e(sigma_v)
        scalar `mu'=e(eu)
        scalar `N'=e(N)
        scalar `conv'=.
        if "`method'"=="qle" scalar `conv'=e(qle_converged)
        foreach key in parameter_valid score_valid estimate_valid qle_optimizer_success qle_optimizer_status qle_optimizer_nfev {
            local value_`key'=e(`key')
        }
        foreach key in invalid_reason score_status raw_sigma_u_text raw_sigma_v_text raw_lambda_text qle_optimizer_message {
            local text_`key' `"`e(`key')'"'
        }
    }
    _estimates unhold `held'
    if `code' exit `code'
    quietly summarize `rawu' if `sample', meanonly
    scalar `um'=r(mean)
    quietly count if `sample' & `rawu'<0
    scalar `neg'=r(N)
    quietly summarize `score' if `sample', meanonly
    scalar `sm'=r(mean)
    preserve
    if "`replace'"!="" {
        foreach v of local destinations {
            capture confirm variable `v', exact
            if !_rc drop `v'
        }
    }
    rename `score' `generate'_score
    rename `rawu' `generate'_raw_u
    label variable `generate'_score "pyStoNED 0.7.3 DDF reference ratio; `method'"
    label variable `generate'_raw_u "Literal public conditional quantity; not geometric DDF target"
    return scalar N=scalar(`N')
    return scalar sigma_u=scalar(`su')
    return scalar sigma_v=scalar(`sv')
    return scalar mu=scalar(`mu')
    return scalar U_mean=scalar(`um')
    return scalar score_mean=scalar(`sm')
    return scalar negative_U_n=scalar(`neg')
    return scalar denominator_min=scalar(`denominator')
    return scalar mean_distance=scalar(`center')
    return scalar converged=scalar(`conv')
    foreach key in parameter_valid score_valid estimate_valid qle_optimizer_success qle_optimizer_status qle_optimizer_nfev {
        return scalar `key'=`value_`key''
    }
    foreach key in invalid_reason score_status raw_sigma_u_text raw_sigma_v_text raw_lambda_text qle_optimizer_message {
        return local `key' `"`text_`key''"'
    }
    return local invalid_policy "`invalid'"
    return local source_status ddf_rows_verified
    return local reference_version "0.7.3"
    return local sign_status "literal_public_distance_residual"
    return local residual_convention "raw_distance_d"
    if "`aligned'"!="" return local residual_convention "counterfactual_centered_minus_d"
    return local implementation pystoned
    return local method "`method'"
    return local generate "`generate'"
    return local cmd stonedddf
    restore, not
    di as text "Reference DDF uses raw distance d; directional production residual is -d. Use diagnose to compare definitions."
    if scalar(`neg')>0 di as text "Warning: literal public conditional quantity contains negative values."
end
