program define _stonedddf_qle_raw, rclass
    version 16
    syntax varname [if] [in], COST(integer) PYTHONExec(string)
    if !inlist(`cost',0,1) exit 198
    local frontier=cond(`cost',"cost","production")
    tempvar score rawu
    tempname held
    _estimates hold `held', restore nullok copy
    capture noisily _stoned_fit `varlist' `if' `in', method(qle) ///
        frontier(`frontier') cet(mult) invalid(raw) pythonexec(`"`pythonexec'"') ///
        generate(`score') inefficiency(`rawu')
    local code=_rc
    if !`code' {
        local su=e(sigma_u)
        local sv=e(sigma_v)
        local conv=e(qle_optimizer_success)
        local lambda=real("`e(raw_lambda_text)'")
    }
    _estimates unhold `held'
    if `code' exit `code'
    return scalar sigma_u=`su'
    return scalar sigma_v=`sv'
    return scalar converged=`conv'
    return scalar lambda=`lambda'
end
