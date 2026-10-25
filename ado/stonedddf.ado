*! version 0.2.1-rawdiag1 06oct2026
program define stonedddf, rclass sortpreserve
    version 16
    syntax , [GENerate(name) TARgets(name) NOISEAdjusted(name) REPLACE ///
        METHOD(string) MAXiter(integer 200) LTOLerance(real 1e-8) ///
        IMPLEMENTation(string) PYTHONExec(string) DIAGnose INVALID(string)]
    if "`implementation'"=="" local implementation pystoned
    local implementation=lower("`implementation'")
    if "`method'"=="" local method mm
    local method=lower("`method'")
    if !inlist("`implementation'","pystoned","legacy") exit 198
    if "`invalid'"=="" local invalid refuse
    local invalid=lower("`invalid'")
    if !inlist("`invalid'","refuse","raw") exit 198
    if "`implementation'"=="legacy" & "`invalid'"=="raw" {
        di as error "invalid(raw) is available only with implementation(pystoned)"
        exit 198
    }
    if !inlist("`method'","mm","qle","psl") exit 198
    if "`implementation'"=="pystoned" & "`method'"=="psl" exit 198
    if "`implementation'"=="pystoned" & "`targets'`noiseadjusted'"!="" {
        di as error "literal public scoring has no geometric targets; use implementation(legacy)"
        exit 198
    }
    if missing(`maxiter',`ltolerance') | `maxiter'<1 | `ltolerance'<=0 | `ltolerance'>.001 exit 198
    if "`generate'"=="" local generate _stonedddf
    local legacy method(`method') generate(`generate') maxiter(`maxiter') ltolerance(`ltolerance') `replace'
    if "`targets'"!="" local legacy `legacy' targets(`targets')
    if "`noiseadjusted'"!="" local legacy `legacy' noiseadjusted(`noiseadjusted')
    if `"`pythonexec'"'=="" local pythonexec python3
    tempname D V pp
    if "`diagnose'"!="" {
        matrix `D'=J(7,12,.)
        matrix colnames `D'=rc sigma_u sigma_v mu U_mean score_mean negative_U_n denominator_min mean_distance converged reference_sign parameter_valid
        matrix rownames `D'=public_mm public_qle aligned_mm aligned_qle legacy_mm legacy_psl legacy_qle
        matrix `V'=J(4,3,.)
        matrix rownames `V'=public_mm public_qle aligned_mm aligned_qle
        matrix colnames `V'=parameter_valid score_valid estimate_valid
        forvalues j=1/7 {
            local meth mm
            if inlist(`j',2,4,7) local meth qle
            if `j'==6 local meth psl
            if `j'<=4 {
                local align
                if inlist(`j',3,4) local align aligned
                capture quietly _stonedddf_public, generate(`pp') method(`meth') ///
                    pythonexec(`"`pythonexec'"') `align' invalid(`invalid')
            }
            else {
                capture quietly _stonedddf_diag_legacy, generate(`pp') method(`meth') ///
                    maxiter(`maxiter') ltolerance(`ltolerance')
            }
            local code=_rc
            matrix `D'[`j',1]=`code'
            matrix `D'[`j',11]=cond(`j'<=2,1,-1)
            if !`code' {
                matrix `D'[`j',12]=1
                if `j'<=4 {
                    matrix `D'[`j',12]=r(parameter_valid)
                    matrix `V'[`j',1]=r(parameter_valid)
                    matrix `V'[`j',2]=r(score_valid)
                    matrix `V'[`j',3]=r(estimate_valid)
                }
                local k=1
                foreach key in sigma_u sigma_v mu U_mean score_mean negative_U_n denominator_min mean_distance converged {
                    local ++k
                    matrix `D'[`j',`k']=r(`key')
                }
            }
            if `code'==498 & inlist(`j',2,4) {
                local align
                if `j'==4 local align aligned
                capture quietly _stonedddf_qle_probe, pythonexec(`"`pythonexec'"') `align'
                if !_rc {
                    matrix `D'[`j',2]=r(sigma_u)
                    matrix `D'[`j',3]=r(sigma_v)
                    matrix `D'[`j',4]=r(sigma_u)*sqrt(2/_pi)
                    matrix `D'[`j',10]=r(converged)
                    matrix `D'[`j',12]=(r(sigma_u)>=0 & r(sigma_v)>0 & !missing(r(sigma_u),r(sigma_v)))
                    matrix `V'[`j',1]=`D'[`j',12]
                    matrix `V'[`j',3]=0
                }
            }
            foreach suffix in score raw_u u v distance {
                capture drop `pp'_`suffix'
            }
        }
    }
    // Diagnostics never select the final estimator or act as a fallback.
    if "`implementation'"=="legacy" {
        _stonedddf_legacy, `legacy'
    }
    else {
        _stonedddf_public, generate(`generate') method(`method') ///
            pythonexec(`"`pythonexec'"') `replace' invalid(`invalid')
    }
    return add
    return local implementation "`implementation'"
    return local invalid_policy "`invalid'"
    if "`diagnose'"!="" {
        return local diagnose_scope "raw public d; counterfactual centered -d; unchanged directional legacy"
        di as text "DDF diagnostics: aligned rows change the residual definition; they are not literal public API results."
        matrix list `D', format(%11.5g)
        return matrix diagnose=`D'
        return matrix diagnose_validity=`V'
        return local diagnose_invalid_policy "`invalid'"
    }
end

program define _stonedddf_diag_legacy, rclass
    version 16
    syntax , GENerate(name) METHod(string) MAXiter(integer) LTOLerance(real)
    _stonedddf_legacy, generate(`generate') method(`method') maxiter(`maxiter') ltolerance(`ltolerance')
    return add
    tempvar sample weighted mapped
    quietly gen byte `sample'=e(sample)
    quietly summarize `generate'_u if `sample', meanonly
    return scalar U_mean=r(mean)
    return scalar negative_U_n=0
    return scalar score_mean=.
    return scalar denominator_min=.
    if "`method'"=="mm" return scalar converged=.
    local coef `e(coefficient_prefix)'
    if "`coef'"!="" {
        quietly gen double `weighted'=0 if `sample'
        local j=0
        local outputs `e(outputs)'
        foreach y of local outputs {
            local ++j
            quietly replace `weighted'=`weighted'+`coef'_gamma`j'*`y' if `sample'
        }
        quietly summarize `weighted' if `sample', meanonly
        local den=r(min)
        return scalar denominator_min=`den'
        quietly count if `sample' & (missing(`weighted') | `weighted'<=0)
        if r(N)==0 {
            quietly gen double `mapped'=1-`generate'_u/`weighted' if `sample'
            quietly summarize `mapped' if `sample', meanonly
            return scalar score_mean=r(mean)
        }
    }
end
