*! version 0.1.0-crs2 01oct2026
*! Joint-input/output CNLS directional distance, VRS/CRS, desirable outputs.
program define cnlsddf_sparse, eclass sortpreserve
    version 16.0
    syntax [if] [in], INputs(varlist numeric min=1) OUTputs(varlist numeric min=1) ///
        GX(string) GY(string) PYthon(string) ///
        [DISTance(name) RESidual(name) COEFficients(name) REPLACE ///
         TOLerance(real 1e-6) NORMTolerance(real 1e-8) MAXiter(integer 200000) ///
         MAXCut(integer 100) PIVot(integer 0) FULL FULLStart(string) TRACE(string) ///
         SOLVer(string) IPMMaxiter(integer 100) CUTBatch(integer 1) INNERCap(real 1e-3) INNERPolicy(string) WORKINGSet(string) RTS(string)]
    // Validate the existing numlist grammar, but do not round explicit literals.
    // A direction containing range notation retains Stata's numlist expansion.
    foreach axis in gx gy {
        capture numlist `"``axis''"', min(1)
        local direction_rc=_rc
        if `direction_rc' exit `direction_rc'
        local expanded "`r(numlist)'"
        mata: _cnlsddfs_direction("`axis'", "expanded", "`axis'_syntax")
    }
    if "`rts'"=="" local rts vrs
    if !inlist("`rts'","vrs","crs") {
        di as error "rts() must be vrs or crs"
        exit 198
    }
    if "`solver'"=="" local solver osqp
    if !inlist("`solver'","osqp","clarabel") {
        di as error "solver() must be osqp or clarabel"
        exit 198
    }
    if "`fullstart'"=="" local fullstart auto
    if !inlist("`fullstart'","auto","cold","interior") | ("`fullstart'"=="interior" & "`full'"=="") {
        di as error "fullstart() must be auto/cold, or interior with full"
        exit 198
    }
    if "`solver'"=="clarabel" & "`fullstart'"!="auto" {
        di as error "fullstart(cold|interior) applies only to solver(osqp)"
        exit 198
    }
    if missing(`cutbatch') | `cutbatch'<1 | missing(`innercap') | `innercap'<=0 | `innercap'>1e-3 {
        di as error "cutbatch() must be positive; innercap() must be in (0, 1e-3]"
        exit 198
    }
    if "`innerpolicy'"=="" local innerpolicy adaptive
    if !inlist("`innerpolicy'","adaptive","monotone") {
        di as error "innerpolicy() must be adaptive or monotone"
        exit 198
    }
    if (`cutbatch'!=1 | `innercap'!=1e-3 | "`innerpolicy'"!="adaptive") & ("`solver'"!="clarabel" | "`full'"!="") {
        di as error "nondefault cutbatch()/innercap()/innerpolicy() require solver(clarabel) without full"
        exit 198
    }
    if "`workingset'"=="" local workingset none
    if !inlist("`workingset'","none","once","bounded") {
        di as error "workingset() must be none, once or bounded"
        exit 198
    }
    if "`workingset'"!="none" & ("`solver'"!="clarabel" | "`full'"!="") {
        di as error "workingset(once|bounded) requires solver(clarabel) without full"
        exit 198
    }
    if "`rts'"=="crs" & ("`solver'"!="clarabel" | "`workingset'"!="none") {
        di as error "rts(crs) requires solver(clarabel) and workingset(none) in this candidate"
        exit 198
    }
    if `"`trace'"'!="" & "`replace'"=="" confirm new file `"`trace'"'
    if !inlist("`c(os)'", "MacOSX", "Unix") {
        di as error "This external-Python bridge currently supports macOS and Unix only"
        exit 198
    }
    marksample touse, novarlist
    markout `touse' `inputs' `outputs'
    local m : word count `inputs'
    local s : word count `outputs'
    local ng : word count `gx'
    local nh : word count `gy'
    local allvars `inputs' `outputs'
    local unique : list uniq allvars
    local nn : word count `unique'
    if `nn' != `m'+`s' {
        di as error "input and output variables must be distinct"
        exit 198
    }
    if `ng'!=`m' | `nh'!=`s' {
        di as error "gx()/gy() must match the input/output dimensions"
        exit 198
    }
    local direction `gx' `gy'
    local positive=0
    foreach g of local direction {
        if missing(`g') | `g'<0 {
            di as error "direction components must be finite and nonnegative"
            exit 198
        }
        if `g'>0 local positive=1
    }
    if !`positive' {
        di as error "the direction cannot be all zero"
        exit 198
    }
    foreach value in tolerance normtolerance {
        if missing(``value'') | ``value''<=0 {
            di as error "`value'() must be finite and positive"
            exit 198
        }
    }
    if missing(`maxiter') | missing(`maxcut') | missing(`ipmmaxiter') | `maxiter'<1 | `maxcut'<1 | `ipmmaxiter'<1 {
        di as error "maxiter(), ipmmaxiter() and maxcut() must be positive integers"
        exit 198
    }
    if missing(`pivot') | `pivot'<0 | `pivot'>`m'+`s' {
        di as error "pivot() must be 0 (automatic) or a direction component index"
        exit 198
    }
    if `pivot'>0 {
        local gp : word `pivot' of `direction'
        if `gp'<=0 {
            di as error "pivot() must select a positive direction component"
            exit 198
        }
    }
    quietly count if `touse'
    local n=r(N)
    if `n'<2 error 2001
    foreach v of local allvars {
        quietly count if `v'<0 & `touse'
        if r(N) {
            di as error "physical inputs and outputs must be nonnegative: `v'"
            exit 459
        }
    }
    if "`distance'"=="" local distance _cnlsddf_distance
    if "`residual'"=="" local residual _cnlsddf_resid
    local destinations `distance' `residual'
    if "`coefficients'"!="" {
        local destinations `destinations' `coefficients'_alpha
        forvalues j=1/`m' {
            local destinations `destinations' `coefficients'_beta`j'
        }
        forvalues j=1/`s' {
            local destinations `destinations' `coefficients'_gamma`j'
        }
    }
    local distinct : list uniq destinations
    local count1 : word count `destinations'
    local count2 : word count `distinct'
    if `count1'!=`count2' {
        di as error "generated variable names must be distinct"
        exit 198
    }
    local overlap : list destinations & allvars
    if "`overlap'"!="" {
        di as error "generated variables cannot replace input/output variables"
        exit 198
    }
    foreach v of local destinations {
        confirm name `v'
        if strlen("`v'")>32 error 198
        if "`replace'"=="" confirm new variable `v'
    }
    quietly findfile cnlsddf_osqp.py
    local backend `"`r(fn)'"'
    confirm file `"`python'"'
    tempfile input bridge backendlog
    // Double-quoted POSIX arguments; reject characters with shell expansion semantics.
    foreach path in python backend input bridge backendlog {
        mata: st_local("unsafe", strofreal(any(strpos(st_local("`path'"), (char(34),char(36),char(96),char(92),char(10),char(13))) :> 0)))
        if `unsafe' {
            di as error "Unsupported shell metacharacter in `path' path"
            exit 198
        }
    }
    mata: _cnlsddfs_export(st_local("input"), "`allvars'", "`touse'")
    local fullarg
    if "`full'"!="" local fullarg --full
    local tracearg
    if `"`trace'"'!="" local tracearg --trace
    quietly shell "`python'" "`backend'" --rts `rts' --input "`input'" --inputs `inputs' --outputs `outputs' --gx `gx' --gy `gy' --output "`bridge'" --tolerance `tolerance' --normtolerance `normtolerance' --maxiter `maxiter' --maxcut `maxcut' --pivot `pivot' `fullarg' --full-start `fullstart' --solver `solver' --ipmmaxiter `ipmmaxiter' --cutbatch `cutbatch' --innercap `innercap' --inner-policy `innerpolicy' --working-set `workingset' `tracearg' >"`backendlog'" 2>&1
    foreach suffix in csv meta history reason json {
        capture confirm file "`bridge'.`suffix'"
        if _rc {
            capture type "`backendlog'"
            di as error "Sparse backend failed; no estimation results were posted"
            foreach cleanup in csv meta history reason json {
                capture erase "`bridge'.`cleanup'"
            }
            exit 498
        }
    }
    local temporary
    forvalues j=1/`count1' {
        tempvar result`j'
        quietly gen double `result`j''=.
        local temporary `temporary' `result`j''
    }
    tempname summary history xdirection ydirection
    capture noisily mata: _cnlsddfs_import(st_local("bridge"),"`temporary'","`touse'",`n',`m'+`s',"`summary'","`history'")
    local backend_rc=_rc
    if !`backend_rc' & `"`trace'"'!="" {
        capture copy "`bridge'.json" `"`trace'"', `replace'
        local backend_rc=_rc
    }
    foreach suffix in csv meta history reason json {
        capture erase "`bridge'.`suffix'"
    }
    if `backend_rc' exit `backend_rc'
    mata: st_matrix("`xdirection'",strtoreal(tokens("`gx'"))); st_matrix("`ydirection'",strtoreal(tokens("`gy'")))
    // Compute row binding before committing any generated destination.
    quietly _cnlsddf_signature `allvars' `temporary', sample(`touse') gx(`xdirection') gy(`ydirection')
    local source_signature "`r(signature)'"
    forvalues j=1/`count1' {
        local dest : word `j' of `destinations'
        local tmp : word `j' of `temporary'
        if "`replace'"!="" {
            capture confirm variable `dest', exact
            if !_rc drop `dest'
        }
        rename `tmp' `dest'
    }
    if "`rts'"=="vrs" label variable `distance' "Centered CNLS-DDF signed distance (sparse QP)"
    else label variable `distance' "CRS CNLS-DDF signed distance (not centered)"
    label variable `residual' "Production residual: negative CNLS-DDF distance"
    ereturn post, esample(`touse')
    local j=0
    foreach name of local sparse_names {
        local ++j
        ereturn scalar `name'=`summary'[1,`j']
    }
    matrix colnames `history'=iteration active solver_status max_violation candidates inner_iters inner_tol
    ereturn matrix cp_history=`history'
    matrix colnames `xdirection'=`inputs'
    matrix colnames `ydirection'=`outputs'
    ereturn matrix gx=`xdirection'
    ereturn matrix gy=`ydirection'
    ereturn scalar resid_sign=-1
    ereturn local stop_reason "`sparse_reason'"
    ereturn local distance_var "`distance'"
    ereturn local resid_var "`residual'"
    ereturn local coefficient_prefix "`coefficients'"
    ereturn local inputs "`inputs'"
    ereturn local outputs "`outputs'"
    ereturn local rts "`rts'"
    ereturn scalar distance_centered=("`rts'"=="vrs")
    ereturn local engine "sparse_`solver'"
    ereturn local qp_solver "`solver'"
    ereturn local working_set "`workingset'"
    ereturn local inner_policy "`innerpolicy'"
    ereturn local estimator "cnlsddf_sparse"
    ereturn local backend_version "0.1.0-crs2"
    ereturn local bridge_version "0.1.0-dir1"
    if "`rts'"=="crs" ereturn local kkt_profile "crs2_reg6_refine"
    ereturn local gx_syntax "`gx_syntax'"
    ereturn local gy_syntax "`gy_syntax'"
    local effective_start cold
    if "`full'"!="" & "`fullstart'"!="cold" local effective_start interior
    if "`solver'"=="clarabel" local effective_start none
    ereturn local full_start "`effective_start'"
    ereturn local initializer "none"
    if e(initializer_seconds)>0 ereturn local initializer "clarabel"
    ereturn local trace_file `"`trace'"'
    // Same statistical model and postestimation contract, different numerical engine.
    if "`rts'"=="vrs" ereturn local source_schema "ddf_rowhash_v1"
    else ereturn local source_schema "ddf_crs_rowhash_v1"
    ereturn local source_signature "`source_signature'"
    ereturn local source_variables "`allvars' `destinations'"
    ereturn local cmd "cnlsddf"
    di as text "Sparse CNLS-DDF (`rts'): N=" as result e(N) as text " RSS=" as result %12.6g e(rss)
    di as text "Convergence=" as result e(converged) as text " (`e(stop_reason)')" ///
        "  Afriat=" as result %9.2e e(max_viol) as text "  normalization=" as result %9.2e e(max_normalization)
    if e(initializer_seconds)>0 di as text "Full-QP initializer (Clarabel): used=" ///
        as result e(initializer_used) as text "  iterations=" as result e(initializer_iterations)
    if e(recovery_attempts)>0 di as text "Working-set recovery: restored=" ///
        as result e(recovery_restored) as text "  attempts=" as result e(recovery_attempts)
    if !e(converged) di as error "Final iterate is not certified; inspect e() diagnostics."
end

version 16.0
mata:
mata set matastrict on
void _cnlsddfs_direction(string scalar name, string scalar expanded,
                         string scalar form)
{
    real rowvector values
    values=strtoreal(tokens(st_local(name)))
    if (cols(values) & !hasmissing(values)) st_local(form,"literal")
    else {
        values=strtoreal(tokens(st_local(expanded)))
        st_local(form,"numlist")
    }
    if (!cols(values) | hasmissing(values)) _error(198,"Finite direction values required")
    st_local(name,invtokens(strtrim(strofreal(values,"%24.17g"))))
}
void _cnlsddfs_export(string scalar path, string scalar vars, string scalar sample)
{
    real scalar f, i
    real matrix z
    z=st_data(.,tokens(vars),sample)
    f=fopen(path,"w")
    fput(f,invtokens(tokens(vars),","))
    for (i=1;i<=rows(z);i++) fput(f,invtokens(strtrim(strofreal(z[i,.],"%24.17g")),","))
    fclose(f)
}
real matrix _cnlsddfs_csv(string scalar path, real scalar header)
{
    real scalar f
    real matrix z
    string matrix line
    f=fopen(path,"r")
    if (header) line=fget(f)
    z=J(0,0,.)
    while ((line=fget(f))!=J(0,0,"")) {
        if (rows(z)==0) z=strtoreal(tokens(subinstr(line,","," ")))
        else z=z\strtoreal(tokens(subinstr(line,","," ")))
    }
    fclose(f)
    if (!rows(z) | hasmissing(z)) _error(498,"Invalid sparse backend numeric output")
    return(z)
}
void _cnlsddfs_import(string scalar stem, string scalar destinations, string scalar sample,
    real scalar n, real scalar k, string scalar summary, string scalar history)
{
    real scalar f
    real matrix z, meta, hist
    string scalar reason, names
    z=_cnlsddfs_csv(stem+".csv",0)
    meta=_cnlsddfs_csv(stem+".meta",1)
    hist=_cnlsddfs_csv(stem+".history",1)
    f=fopen(stem+".meta","r"); names=invtokens(tokens(subinstr(fget(f),","," "))); fclose(f)
    f=fopen(stem+".reason","r"); reason=strtrim(fget(f)); fclose(f)
    if (rows(z)!=n | cols(z)!=k+3 | rows(meta)!=1 | cols(meta)!=cols(tokens(names)) | cols(hist)!=7) _error(498,"Sparse backend dimension mismatch")
    if (!any(reason:==("certified","maxcut","inner_maxiter","solver_failure"))) _error(498,"Invalid sparse backend status")
    if (max(abs(z[,1]+z[,2]))!=0) _error(498,"Sparse backend residual sign mismatch")
    st_matrix(summary,meta)
    st_matrix(history,hist)
    st_local("sparse_names",names)
    st_local("sparse_reason",reason)
    st_store(.,tokens(destinations),sample,z[,1..cols(tokens(destinations))])
}
end
