*! version 0.4.3-integrated1 06oct2026
program define stoned, eclass sortpreserve
    version 16
    syntax varname(max=1 numeric) [if] [in], ///
        [FRONTier(string) METHod(string) GENerate(name) INEFFiciency(name) ///
        IMPLEMENTation(string) CET(string) OUTcome(varname numeric) ///
        PYTHONExec(string) VERBose DIAGnose INVALID(string)]
    local residual `varlist'
    if "`diagnose'" == "" {
        _stoned_fit `0'
        exit
    }
    if "`invalid'"=="" local invalid refuse
    local invalid=lower("`invalid'")
    if !inlist("`invalid'","refuse","raw") exit 198
    if lower("`implementation'")=="legacy" & "`invalid'"=="raw" exit 198
    marksample touse
    if "`implementation'" == "" local implementation "pystoned"
    local implementation=lower("`implementation'")
    if "`frontier'" == "" local frontier "production"
    local frontier=lower("`frontier'")
    if "`cet'" == "" local cet "addi"
    local cet=lower("`cet'")
    if !inlist("`frontier'","production","cost") | ///
       !inlist("`implementation'","pystoned","legacy") | ///
       !inlist("`cet'","addi","mult") exit 198
    if "`cet'" == "addi" & "`outcome'" == "" {
        di as error "diagnose with ADDI comparison requires outcome(y); cet(mult) compares MULT scores."
        exit 198
    }
    if "`outcome'" != "" {
        quietly count if `touse' & (missing(`outcome') | `outcome'==0)
        if r(N) exit 498
    }
    if "`cet'" == "mult" & "`outcome'" != "" exit 198
    local output_score "`generate'"
    if "`output_score'"=="" local output_score "_stoned_te"
    local output_names `output_score' `inefficiency'
    local protected `residual' `outcome'
    local collisions : list output_names & protected
    if "`collisions'"!="" exit 198
    local is_cost=("`frontier'"=="cost")
    local sign=cond(`is_cost',1,-1)
    if "`pythonexec'" == "" local pythonexec "python3"

    // Diagnostics use temporary outputs. Restore prior e() before the selected fit.
    tempname held D V K rawsu rawsv rawconv rawlam
    _estimates hold `held', restore nullok copy
    matrix `D'=J(4,17,.)
    matrix colnames `D'=rc sigma_u sigma_v lambda valid converged nll_mean moment_gap U_mean BC_mean mapped_mean literal_mean max_U_gap delta_zero literal_U_min negative_U_n flat_zero
    matrix rownames `D'=public_mm public_qle legacy_mm legacy_psl
    matrix `V'=J(2,3,.)
    matrix rownames `V'=public_mm public_qle
    matrix colnames `V'=parameter_valid score_valid estimate_valid
    tempvar dte du
    local psl_effective "failed"
    local psl_selection "failed"
    forvalues j=1/4 {
        local impl "legacy"
        local meth "mm"
        local extra ""
        if `j'==1 {
            local impl "pystoned"
            local extra "cet(`cet') invalid(`invalid')"
            if "`outcome'"!="" local extra "`extra' outcome(`outcome')"
        }
        if `j'==2 & "`invalid'"=="raw" {
            local impl pystoned
            local meth qle
            local extra `"cet(`cet') invalid(raw) pythonexec("`pythonexec'")"'
            if "`outcome'"!="" local extra "`extra' outcome(`outcome')"
        }
        if `j'==4 local meth "psl"
        if `j'!=2 | "`invalid'"=="raw" {
            capture quietly _stoned_fit `residual' if `touse', ///
                implementation(`impl') method(`meth') frontier(`frontier') ///
                generate(`dte') inefficiency(`du') `extra'
            local code=_rc
            matrix `D'[`j',1]=`code'
            if !`code' {
                matrix `D'[`j',2]=e(sigma_u)
                matrix `D'[`j',3]=e(sigma_v)
                matrix `D'[`j',4]=e(lambda)
                matrix `D'[`j',5]=1
                if `j'<=2 {
                    matrix `D'[`j',5]=e(parameter_valid)
                    matrix `V'[`j',1]=e(parameter_valid)
                    matrix `V'[`j',2]=e(score_valid)
                    matrix `V'[`j',3]=e(estimate_valid)
                    if `j'==2 matrix `D'[`j',6]=e(qle_optimizer_success)
                }
                if `j'==4 {
                    matrix `D'[`j',6]=e(psl_converged)
                    local psl_effective "`e(effective_method)'"
                    local psl_selection "`e(psl_selection)'"
                }
            }
        }
        else {
            // Raw public QLE permits observing invalid parameters and optimizer success.
            // The refuse-policy probe does not supply selected score outputs.
            capture quietly _stoned_qle_raw `residual' if `touse', ///
                cost(`is_cost') pythonexec(`"`pythonexec'"')
            local code=_rc
            matrix `D'[`j',1]=`code'
            if !`code' {
                matrix `D'[`j',2]=r(sigma_u)
                matrix `D'[`j',3]=r(sigma_v)
                matrix `D'[`j',4]=r(lambda)
                matrix `D'[`j',5]=(r(sigma_u)>=0 & r(sigma_v)>0 & !missing(r(sigma_u),r(sigma_v)))
                matrix `D'[`j',6]=r(converged)
                if !`D'[`j',5] | !`D'[`j',6] matrix `D'[`j',1]=498
                matrix `V'[2,1]=r(parameter_valid)
                matrix `V'[2,2]=.
                matrix `V'[2,3]=cond(`D'[2,1]==498,0,.)
            }
        }
        if `D'[`j',5]==1 {
            // All methods are compared using the same score definition and residuals.
            mata: st_matrix("`K'",_stoned_diag_metrics("`residual'","`touse'", ///
                "`outcome'","`cet'",st_matrix("`D'")[`j',2], ///
                st_matrix("`D'")[`j',3],`sign'))
            forvalues c=1/10 {
                matrix `D'[`j',`c'+6]=`K'[1,`c']
            }
            if `j'==2 matrix `D'[`j',17]=(abs(`D'[`j',14])<=1e-8 & abs(`D'[`j',4])<=.1)
        }
    }
    quietly summarize `residual' if `touse', meanonly
    local diag_mean=r(mean)
    local diag_N=r(N)
    _estimates unhold `held'
    di as text _n "StoNED diagnostic comparison: fixed JSS pyStoNED 0.7.3 / retained RC3"
    di as text "Same residual sample; lower common NLL is a fit diagnostic, not accuracy or model selection."
    matrix list `D', format(%11.5g) noheader
    di as text "MOM: public uses divisor N; legacy uses N-1 with corresponding third-moment scaling."
    di as text "QLE constrains total variance by raw mean(residual^2); PSL estimates two scales."
    di as text "BC_mean and mapped_mean have different definitions; mapped_mean uses stable U with the same `cet' mapping."
    di as text "literal_mean uses the published denominator; max_U_gap isolates its difference from stable U at SAME parameters."
    if `D'[2,5]==0 {
        di as text "QLE PARAMETER INVALID: optimizer success does not make negative sigma_u admissible."
    }
    if `D'[2,17]==1 {
        di as text "QLE NEAR-ZERO FLATNESS: near-boundary parameter changes can have negligible objective change."
        di as text "This numerical screen is not a test of zero inefficiency or a confidence interval."
    }
    if `D'[2,1]!=0 & missing(`D'[2,5]) {
        di as text "PUBLIC QLE UNAVAILABLE: verify pythonexec(), pinned runtime and source hashes. Other rows remain visible."
    }
    forvalues j=1/4 {
        if `D'[`j',16]>0 & `D'[`j',16]<. {
            di as text "REFERENCE TAIL INVALID: row `j' yields negative conditional U under the literal public formula."
        }
    }
    di as text "Requested legacy PSL: effective method=`psl_effective'; selection=`psl_selection'."
    di as text "Input residual mean=" as result %12.5g `diag_mean' as text "; N=" as result `diag_N'
    di as text "Residual levels are not recentered. RC=498 alone does not identify a software defect."
    di as text "No method is automatically selected. Failed rows are retained; they do not change the requested fit."

    local opts "implementation(`implementation') frontier(`frontier') invalid(`invalid')"
    if "`method'"!="" local opts "`opts' method(`method')"
    if "`generate'"!="" local opts "`opts' generate(`generate')"
    if "`inefficiency'"!="" local opts "`opts' inefficiency(`inefficiency')"
    if "`implementation'"=="pystoned" {
        local opts "`opts' cet(`cet')"
        if "`outcome'"!="" local opts "`opts' outcome(`outcome')"
        if lower("`method'")=="qle" local opts `"`opts' pythonexec("`pythonexec'")"'
    }
    // Leave exact selected-fit errors intact. No diagnostics-as-success substitution.
    _stoned_fit `residual' `if' `in', `opts' `verbose'
    ereturn matrix diagnose=`D'
    ereturn matrix diagnose_validity=`V'
    ereturn local diagnose_invalid_policy "`invalid'"
    ereturn local diagnose_reference "JSS-pyStoNED-0.7.3"
    ereturn local diagnose_score "`cet'"
    ereturn local diagnose_psl_effective "`psl_effective'"
    ereturn local diagnose_psl_selection "`psl_selection'"
    ereturn scalar diagnose_resid_mean=`diag_mean'
    ereturn scalar diagnose_N=`diag_N'
    ereturn scalar diagnose_flat_tol=1e-8
    ereturn local diagnose_policy "comparison_only_no_selection"
end

program define _stoned_qle_raw, rclass
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
    if !`code' {
        local pv=e(parameter_valid)
        local svflag=e(score_valid)
        local ev=e(estimate_valid)
    }
    _estimates unhold `held'
    if `code' exit `code'
    return scalar parameter_valid=`pv'
    return scalar score_valid=`svflag'
    return scalar estimate_valid=`ev'
    return scalar sigma_u=`su'
    return scalar sigma_v=`sv'
    return scalar converged=`conv'
    return scalar lambda=`lambda'
end

mata:
real scalar _stoned_diag_tail_excess(real scalar t)
{
    real scalar k, tail
    tail = 0
    for (k=100; k>=2; k--) tail = k/(t+tail)
    return(1/(t+tail))
}

// Conditional N(mu, ss^2) truncated below at zero: (E(exp(-u)), E(u)).
// Negative tails use Mills ratios; the mixed-tail identity avoids ss^2 overflow.
real rowvector _stoned_diag_conditional(real scalar mu, real scalar ss)
{
    real scalar a, b, t, q, logte, u
    if (hasmissing((mu,ss)) | ss<=0) return((.,.))
    a = mu/ss
    if (missing(a)) return((.,.))
    b = a-ss
    if (a < -8) {
        t = -a
        q = _stoned_diag_tail_excess(t)
        u = ss*q
        logte = -ln(1+(ss+_stoned_diag_tail_excess(t+ss)-q)/(t+q))
    }
    else {
        u = mu + ss*exp(lnnormalden(a)-lnnormal(a))
        if (b < -8) {
            t = -b
            logte = -.5*a^2 -.5*ln(2*pi()) - ///
                ln(t+_stoned_diag_tail_excess(t)) - lnnormal(a)
        }
        else logte = -mu + ss^2/2 + lnnormal(b)-lnnormal(a)
    }
    if (hasmissing((logte,u)) | logte>1e-12 | u<0) return((.,.))
    return((exp(min((logte,0))),u))
}


real rowvector _stoned_diag_metrics(string scalar rv,string scalar tv,
    string scalar yv,string scalar cet,real scalar su,real scalar sv,real scalar sign)
{
    real colvector eps,z,a,u,bc,pu,ps,ls,y
    real rowvector pair
    real scalar sig,ss,i,nll,zero,raw2,mgap
    eps=st_data(.,rv,tv)
    sig=sqrt(su^2+sv^2)
    ss=su*sv/sig
    z=(eps:+sign*su*sqrt(2/pi())):/sig
    raw2=mean(eps:^2)
    nll=ln(sig)+mean(.5:*z:^2:-lnnormal(sign*su/sv:*z))
    zero=.5*ln(raw2)+.5+ln(2)
    mgap=(sv^2+(1-2/pi())*su^2-raw2)/raw2
    u=J(rows(eps),1,0);bc=J(rows(eps),1,1)
    if (su>0) {
        for(i=1;i<=rows(eps);i++) {
            pair=_stoned_diag_conditional(sign*z[i]*su^2/sig,ss)
            bc[i]=pair[1];u[i]=pair[2]
        }
    }
    a=z:*su/sv
    pu=ss:*(normalden(a):/(1:-normal(-sign:*a):+1e-6):+sign:*a)
    if (cet=="addi") {
        y=st_data(.,yv,tv)
        ps=(y:+sign:*u):/y
        ls=(y:+sign:*pu):/y
    }
    else {
        ps=exp(sign:*u)
        ls=exp(sign:*pu)
    }
    if (hasmissing((u,bc,pu,ps,ls)) | hasmissing((nll,zero,mgap))) return(J(1,10,.))
    return((nll,mgap,mean(u),mean(bc),mean(ps),mean(ls),max(abs(u-pu)),nll-zero,min(pu),sum(pu:<0)))
}
end
