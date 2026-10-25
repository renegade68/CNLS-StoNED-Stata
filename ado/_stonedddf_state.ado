program define _stonedddf_state, rclass
    version 16
    syntax , SAMPLE(name)
    if "`e(cmd)'"!="cnlsddf" | "`e(rts)'"!="vrs" {
        di as error "stonedddf requires the current cnlsddf VRS estimation results"
        exit 301
    }
    if e(converged)!=1 {
        di as error "cnlsddf did not certify convergence; refit before decomposition"
        exit 498
    }
    if "`generate'"=="" local generate _stonedddf
    local inputs `e(inputs)'
    local outputs `e(outputs)'
    local distance `e(distance_var)'
    local residual `e(resid_var)'
    local m=e(m)
    local s=e(s)
    local touse `sample'
    quietly gen byte `touse'=e(sample)
    quietly count if `touse'
    if r(N)!=e(N) {
        di as error "estimation sample changed; refit cnlsddf"
        exit 459
    }
    if r(N)<3 error 2001
    local protected `inputs' `outputs' `distance' `residual'
    foreach v of local protected {
        capture confirm numeric variable `v', exact
        if _rc exit 459
        quietly count if missing(`v') & `touse'
        if r(N) {
            di as error "missing stage-one data in the estimation sample: `v'"
            exit 459
        }
    }
    foreach v of varlist `inputs' `outputs' {
        quietly count if `v'<0 & `touse'
        if r(N) error 459
    }
    local coef `e(coefficient_prefix)'
    if "`coef'"!="" {
        local protected `protected' `coef'_alpha
        forvalues j=1/`m' {
            local protected `protected' `coef'_beta`j'
        }
        forvalues j=1/`s' {
            local protected `protected' `coef'_gamma`j'
        }
    }
    // Old fits lack a joint-row binding: refit explicitly with this candidate.
    if "`e(source_schema)'"!="ddf_rowhash_v1" | "`e(source_signature)'"=="" | ///
        "`e(source_variables)'"!="`protected'" {
        di as error "stage-one row provenance missing or incompatible; refit with the DDF safety candidate"
        exit 459
    }
    foreach v of local protected {
        capture confirm numeric variable `v', exact
        if _rc {
            di as error "stage-one source variable missing or nonnumeric: `v'; refit cnlsddf"
            exit 459
        }
    }
    tempname source_gx source_gy
    matrix `source_gx'=e(gx)
    matrix `source_gy'=e(gy)
    quietly _cnlsddf_signature `protected', sample(`touse') gx(`source_gx') gy(`source_gy')
    if "`r(signature)'"!="`e(source_signature)'" {
        di as error "stage-one rows or direction changed; refit cnlsddf"
        exit 459
    }
    quietly count if `touse' & abs(`residual'+`distance')>1e-10*max(1,abs(`distance'))
    if r(N) exit 459
    return local protected "`protected'"
    return local distance "`distance'"
    return local outputs "`outputs'"
    return local coefficient "`coef'"
end
