program define _stonedddf_qle_probe, rclass
    version 16
    syntax , PYTHONExec(string) [ALIGNED]
    tempvar sample eps weighted
    quietly _stonedddf_state, sample(`sample')
    local d `r(distance)'
    local outputs `r(outputs)'
    local coef `r(coefficient)'
    if "`coef'"=="" exit 459
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
    _stonedddf_qle_raw `eps' if `sample', cost(0) pythonexec(`"`pythonexec'"')
    return add
end
