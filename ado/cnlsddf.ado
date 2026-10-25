*! version 0.1.0-safety1 01oct2026
*! Joint-input/output CNLS directional distance, VRS, desirable outputs.
program define cnlsddf, eclass sortpreserve
    version 16.0
    syntax [if] [in], INputs(varlist numeric min=1) OUTputs(varlist numeric min=1) ///
        GX(numlist min=1) GY(numlist min=1) ///
        [DISTance(name) RESidual(name) COEFficients(name) REPLACE ///
         TOLerance(real 1e-6) NORMTolerance(real 1e-8) MAXiter(integer 20000) ///
         MAXCut(integer 100) RHO(real 1) PIVot(integer 0) FULL COLDStart CPState(string) INNERPolicy(string) VERBose]
    if "`innerpolicy'"=="" local innerpolicy adaptive
    local innerpolicy=lower("`innerpolicy'")
    if !inlist("`innerpolicy'","adaptive","strict") {
        di as error "innerpolicy() must be adaptive or strict"
        exit 198
    }
    if "`cpstate'"=="" local cpstate cold
    local cpstate=lower("`cpstate'")
    if !inlist("`cpstate'","cold","mapped") | ("`coldstart'"!="" & "`cpstate'"=="mapped") {
        di as error "cpstate() must be cold or mapped; coldstart cannot be combined with mapped"
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
    foreach g of numlist `direction' {
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
    foreach value in tolerance normtolerance rho {
        if missing(``value'') | ``value''<=0 {
            di as error "`value'() must be finite and positive"
            exit 198
        }
    }
    if missing(`maxiter') | missing(`maxcut') | `maxiter'<1 | `maxcut'<1 {
        di as error "maxiter() and maxcut() must be positive integers"
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
    // Stage results in temporary variables; invalid calls never overwrite targets.
    local temporary
    forvalues j=1/`count1' {
        tempvar result`j'
        quietly gen double `result`j''=.
        local temporary `temporary' `result`j''
    }
    local dtmp : word 1 of `temporary'
    local rtmp : word 2 of `temporary'
    local coeftmp : list temporary - dtmp
    local coeftmp : list coeftmp - rtmp
    mata: _cnlsddf_estimate("`inputs'","`outputs'","`touse'","`gx'","`gy'", ///
        "`dtmp'","`rtmp'","`coeftmp'",`tolerance',`normtolerance',`maxiter', ///
        `maxcut',`rho',`pivot',("`full'"!=""),("`verbose'"!=""),("`cpstate'"=="cold"),("`innerpolicy'"=="adaptive"))
    tempname summary inner history xdirection ydirection
    matrix `summary'=r(ddf_summary)
    matrix `inner'=r(ddf_inner)
    matrix `history'=r(ddf_history)
    matrix `xdirection'=r(ddf_gx)
    matrix `ydirection'=r(ddf_gy)
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
    label variable `distance' "Centered CNLS-DDF signed distance"
    label variable `residual' "Production residual: negative CNLS-DDF distance"
    ereturn post, esample(`touse')
    local names converged cp_iters n_active max_viol max_mono max_mono_scaled max_normalization ///
        rss total_inner_iters inner_solves warm_starts factorizations chol_solves lu_fallbacks ///
        rho_updates final_rho pivot parameters constraint_bytes dense_constraint_bytes mean_distance expanded_starts
    local j=0
    foreach name of local names {
        local ++j
        ereturn scalar `name'=`summary'[1,`j']
    }
    local names inner_status inner_iter inner_primal inner_dual inner_eps_primal inner_eps_dual ///
        inner_stationarity inner_eps_stationarity inner_complementarity inner_eps_complementarity inner_tolerance
    local j=0
    foreach name of local names {
        local ++j
        ereturn scalar `name'=`inner'[1,`j']
    }
    matrix colnames `history'=iteration active inner_status max_violation candidates inner_iters inner_tol
    ereturn matrix cp_history=`history'
    matrix colnames `xdirection'=`inputs'
    matrix colnames `ydirection'=`outputs'
    ereturn matrix gx=`xdirection'
    ereturn matrix gy=`ydirection'
    ereturn scalar N=`n'
    ereturn scalar m=`m'
    ereturn scalar s=`s'
    ereturn scalar tolerance=`tolerance'
    ereturn scalar normtolerance=`normtolerance'
    ereturn scalar resid_sign=-1
    ereturn scalar normal_matrix_bytes=8*e(parameters)^2
    ereturn local stop_reason "`ddf_reason'"
    ereturn local distance_var "`distance'"
    ereturn local resid_var "`residual'"
    ereturn local coefficient_prefix "`coefficients'"
    ereturn local inputs "`inputs'"
    ereturn local outputs "`outputs'"
    ereturn local rts "vrs"
    ereturn local inner_policy "`innerpolicy'"
    ereturn local expansion "`cpstate'"
    ereturn local engine "structured_reduced_admm"
    ereturn local source_schema "ddf_rowhash_v1"
    ereturn local source_signature "`source_signature'"
    ereturn local source_variables "`allvars' `destinations'"
    ereturn local cmd "cnlsddf"
    di as text "CNLS-DDF (VRS): N=" as result e(N) as text " inputs=" as result e(m) ///
        as text " outputs=" as result e(s) as text " RSS=" as result %12.6g e(rss)
    di as text "Convergence=" as result e(converged) as text " (`e(stop_reason)')" ///
        "  Afriat=" as result %9.2e e(max_viol) as text "  normalization=" as result %9.2e e(max_normalization)
    if !e(converged) di as error "Final iterate is not certified; inspect e() diagnostics."
end

// Standalone fork of the dev5 structured ADMM kernel, extended for G*t <= upper.
// Keep convnls unchanged; validate this kernel against dense operators and R QPs.
version 16.0
mata:
mata set matastrict on
struct _cnlsddf_operator {
    real scalar dense, n, p, q, sign, nvars, ncons, raw_rho
    real matrix G, Gt, Z, pairs, zh, zhi, pi, ph
    real rowvector bound
    real colvector upper
    real colvector oi, oh, ui, uh, mono
}

// Continuation belongs to one fixed QP and its returned primal iterate.
// The caller must invalidate it whenever Q, c, or the ordered constraints change.
struct _cnlsddf_admm_state {
    real scalar valid, reused, nvars, ncons, rho, iteration, pending_rho, raw_rho
    real scalar alpha, prev_primal, diverge_count, primal, dual, eps_primal, eps_dual
    real colvector s, u
}

// Balance distance from the existing stopping thresholds, not raw residual units.
// These ratios guide the penalty only; all four certification checks stay intact.
real scalar _cnlsddf_next_rho(real scalar rho, real scalar primal, real scalar dual,
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
real scalar _cnlsddf_penalty(real scalar rho, real scalar primal, real scalar dual,
                           real scalar eps_primal, real scalar eps_dual,
                           real scalar raw_rho)
{
    if (!raw_rho) return(_cnlsddf_next_rho(rho,primal,dual,eps_primal,eps_dual))
    if (primal > 5*dual & dual > 0) return(rho*2)
    if (dual > 5*primal & primal > 0) return(max((rho/2,0.1)))
    return(rho)
}

struct _cnlsddf_operator scalar _cnlsddf_pair_operator(real matrix X,
    real matrix pairs, real scalar q, real scalar is_crs,
    real scalar is_concave, real scalar do_mono, real rowvector bound)
{
    struct _cnlsddf_operator scalar op
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
    op.bound=bound
    op.ncons=a+rows(op.mono)+op.n
    op.upper=(J(a+rows(op.mono),1,0) \ J(op.n,1,1))
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

real colvector _cnlsddf_op_times(struct _cnlsddf_operator scalar op,
                                real colvector d)
{
    real matrix D
    real colvector out
    real scalar a
    if (op.dense) return(op.G*d)
    if (hasmissing(d)) return(J(op.ncons,1,.))
    a=rows(op.pairs)
    D=colshape(d,op.p)
    out=J(op.ncons,1,0)
    if (a) {
        D=colshape(d[|1\op.n*op.p|],op.p)
        out[|1\a|]=op.sign*rowsum(op.zh:*(D[op.pairs[,1],.]-D[op.pairs[,2],.]),1)
    }
    if (rows(op.mono)) out[|a+1\a+rows(op.mono)|]=-d[op.mono]
    out[|a+rows(op.mono)+1\op.ncons|]=D*op.bound'
    return(out)
}

real colvector _cnlsddf_op_transpose(struct _cnlsddf_operator scalar op,
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
    D=op.sign*D + v[|a+rows(op.mono)+1\op.ncons|]*op.bound
    out=(vec(D') \ J(op.q,1,0))
    if (rows(op.mono)) out[op.mono]=out[op.mono]-v[|a+1\a+rows(op.mono)|]
    return(out)
}

real matrix _cnlsddf_op_gram(struct _cnlsddf_operator scalar op)
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
    for (k=1;k<=op.n;k++) {
        ii=(k-1)*op.p:+(1..op.p)
        out[ii,ii]=out[ii,ii]+op.bound'*op.bound
    }
    return(out)
}

real scalar _cnlsddf_op_bytes(struct _cnlsddf_operator scalar op)
{
    // Matrix payload only, excluding allocator overhead and temporary workspaces.
    return(8*(length(op.Z)+length(op.pairs)+length(op.zh)+length(op.zhi)+
              length(op.pi)+length(op.ph)+length(op.oi)+length(op.oh)+
              length(op.ui)+length(op.uh)+length(op.mono)+length(op.bound)+length(op.upper)))
}

// Append-only constraints, fixed objective and returned primal iterate.
// Preserve multipliers for unchanged rows; initialize new rows with zero dual.
// Caller owns Q/c identity. Reject any change in old constraint order or scaling.
real scalar _cnlsddf_expand_state(struct _cnlsddf_operator scalar old,
    struct _cnlsddf_operator scalar next, real colvector t,
    struct _cnlsddf_admm_state scalar state)
{
    real scalar a,b
    real colvector map,slack,dual
    if (state.valid!=1) return(0)
    state.valid=0
    if (old.dense | next.dense) return(0)
    if (old.n!=next.n | old.p!=next.p | old.q!=next.q |
        old.nvars!=next.nvars | old.sign!=next.sign | old.raw_rho!=next.raw_rho) return(0)
    a=rows(old.pairs); b=rows(next.pairs)
    if (b<=a | rows(old.mono)!=rows(next.mono)) return(0)
    if (max(abs(old.Z-next.Z))!=0 | max(abs(old.bound-next.bound))!=0) return(0)
    if (a) if (max(abs(old.pairs-next.pairs[|1,1\a,2|]))!=0) return(0)
    if (rows(old.mono)) if (max(abs(old.mono-next.mono))!=0) return(0)
    if (state.ncons!=old.ncons | state.nvars!=old.nvars | rows(t)!=old.nvars) return(0)
    if (rows(state.s)!=old.ncons | rows(state.u)!=old.ncons) return(0)
    if (hasmissing(state.s) | hasmissing(state.u) | hasmissing(t) |
        missing(state.rho) | state.rho<=0) return(0)
    map=((1::a) \ ((a+1::old.ncons):+(b-a)))
    if (!a) map=(1::old.ncons):+b
    if (max(abs(old.upper-next.upper[map]))!=0) return(0)
    slack=_cnlsddf_op_times(next,t)
    if (hasmissing(slack) | hasmissing(next.upper)) return(0)
    slack=rowmin((slack,next.upper))
    dual=J(next.ncons,1,0)
    slack[map]=state.s; dual[map]=state.u
    state.s=slack; state.u=dual; state.ncons=next.ncons
    // New problem: do not adapt rho using residuals from the previous QP.
    state.iteration=0; state.pending_rho=0
    state.prev_primal=.; state.diverge_count=0
    state.valid=1
    return(1)
}

real colvector _cnlsddf_admm_resume(real matrix ZtZ, real colvector Zty,
                                 struct _cnlsddf_operator scalar op, real colvector delta0,
                                 real scalar rho, real scalar tol,
                                 real scalar maxiter, real scalar verbose,
                                 real rowvector info, struct _cnlsddf_admm_state scalar state)
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
    Gdelta = _cnlsddf_op_times(op, delta)
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
        if (state.pending_rho) {
            new_rho=_cnlsddf_penalty(rho,state.primal,state.dual,
                                      state.eps_primal,state.eps_dual,raw_rho)
            if (new_rho != rho) {
                u=u:*(rho/new_rho)
                rho=new_rho
                info[15]=info[15]+1
            }
        }
    }
    else {
        s = rowmin((Gdelta, op.upper))
        u = J(n_con, 1, 0)
        offset=0
        alpha_relax=1.6
        prev_r_pri=.
        diverge_count=0
    }
    GtG = _cnlsddf_op_gram(op)
    reg = 1e-8
    H = ZtZ + rho * GtG + reg * I(n_vars)
    factor_dirty = 1
    factor_ok = 0
    for (k = 1; k <= maxiter; k++) {
        iteration=offset+k
        delta_prev = delta
        // Proximal stabilization: matched RHS term preserves the unregularized QP.
        rhs = Zty + rho * _cnlsddf_op_transpose(op, s - u) + reg * delta_prev
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
        Gdelta = _cnlsddf_op_times(op, delta)
        Gdelta_hat = alpha_relax * Gdelta + (1-alpha_relax)*s
        s_old = s
        s = rowmin((Gdelta_hat + u, op.upper))
        u = u + Gdelta_hat - s

        // Always measure the last iterate, even when maxiter is not a multiple of 25.
        regular_check=(mod(iteration,25)==0 | iteration<=5)
        if (regular_check | k == maxiter) {
            dual_term = rho * _cnlsddf_op_transpose(op, u)
            grad = ZtZ * delta
            r_pri = max(abs(Gdelta - s))
            r_dual = rho * max(abs(_cnlsddf_op_transpose(op, s - s_old)))
            eps_pri = tol + tol * max((max(abs(Gdelta)), max(abs(s))))
            eps_dual = tol + tol * max(abs(dual_term))
            stationarity = max(abs(grad - Zty + dual_term))
            eps_stat = tol + tol * max((max(abs(grad)), max(abs(Zty)), max(abs(dual_term))))
            comp = max(abs((s-op.upper) :* (rho*u)))
            eps_comp = tol + tol * max(abs(s-op.upper)) * max(abs(rho*u))
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
                s = rowmin((Gdelta, op.upper))
                diverge_count = 0
                prev_r_pri = .
            }
        }
        if (mod(iteration,100) == 0 & k < maxiter) {
            new_rho = _cnlsddf_penalty(rho,r_pri,r_dual,eps_pri,eps_dual,raw_rho)
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
    state.pending_rho=(info[1]==0 & mod(state.iteration,100)==0)
    return(delta)
}

// Full separation in original physical units; pairs are (alternative, own).
// Keep at most one NEW worst violation per own observation; scan every pair.
real matrix _cnlsddf_scan(real matrix W, real matrix coef, real matrix added,
                          real scalar tol, real scalar violation)
{
    real scalar i,h,n,besth,best,v
    real matrix candidates
    real colvector own
    n=rows(W)
    own=rowsum(W:*coef,1)
    candidates=J(0,2,0)
    violation=0
    if (hasmissing(coef) | hasmissing(own)) {
        violation=.
        return(candidates)
    }
    for (i=1;i<=n;i++) {
        best=tol; besth=0
        for (h=1;h<=n;h++) if (h!=i) {
            v=own[i]-W[i,.]*coef[h,.]'
            if (missing(v)) {
                violation=.
                return(candidates)
            }
            violation=max((violation,v))
            if (!added[h,i] & v>best) {
                best=v
                besth=h
            }
        }
        if (besth) candidates=candidates\(besth,i)
    }
    return(candidates)
}

void _cnlsddf_estimate(string scalar xs, string scalar ys, string scalar use,
    string scalar gxstr, string scalar gystr, string scalar dvar, string scalar rvar,
    string scalar coefvars, real scalar tol, real scalar ntol, real scalar maxiter,
    real scalar maxcut, real scalar rho, real scalar pivot, real scalar full,
    real scalar verbose, real scalar coldstart, real scalar adaptive)
{
    real matrix X,Y,V,W,B,Q,T,coeff,slopes,pairs,added,newpairs,hist
    real rowvector scales,g,a,free,idx,info,stats,summary
    real colvector q,c,t,d,obs
    real scalar n,m,s,K,P,i,j,h,cp,conv,viol,mono,monos,normerr,solves,total,warm,inner_tol,rss,expanded
    string scalar reason
    struct _cnlsddf_operator scalar op,nextop
    struct _cnlsddf_admm_state scalar state
    X=st_data(.,xs,use); Y=st_data(.,ys,use)
    n=rows(X); m=cols(X); s=cols(Y); K=m+s; P=n*K
    g=(strtoreal(tokens(gxstr)),strtoreal(tokens(gystr)))
    V=(X,-Y); W=(J(n,1,1),V)
    scales=colmax(abs(V)); for (j=1;j<=K;j++) if (scales[j]==0) scales[j]=1
    V=V:/scales; a=g:/scales
    if (!pivot) {
        idx=selectindex(a:==max(a))
        pivot=idx[1]
    }
    free=select((1..K),(1..K):!=pivot)
    B=(J(n,1,1),V[,free]-V[,pivot]*(a[free]/a[pivot]))
    q=V[,pivot]/a[pivot]
    if (hasmissing(B) | hasmissing(q) | hasmissing(a)) {
        errprintf("direction/quantity scaling overflow; rescale quantities or direction\n")
        _error(498)
    }
    Q=J(P,P,0); c=J(P,1,0)
    for (i=1;i<=n;i++) {
        idx=(i-1)*K:+(1..K)
        Q[idx,idx]=B[i,.]'*B[i,.]
        c[idx]=-B[i,.]'*q[i]
    }
    if (hasmissing(Q) | hasmissing(c)) _error(498)
    T=J(n,K,0); T[,1]=-q
    t=vec(T')
    pairs=J(0,2,0); added=I(n)
    // Adjacent bidirectional pairs are only a starting set, not certification.
    for (i=1;i<=n;i++) for (h=1;h<=n;h++) if (i!=h) {
        if (full | abs(i-h)==1) {
            pairs=pairs\(h,i)
            added[h,i]=1
        }
    }
    state.valid=0; conv=0; solves=total=warm=expanded=0; stats=J(1,4,0)
    hist=J(0,7,0); reason="maxcut"
    inner_tol=tol*.05
    op=_cnlsddf_pair_operator(B[,2..K],pairs,0,0,1,1,(0,a[free]))
    for (cp=1;cp<=maxcut;cp++) {
        t=_cnlsddf_admm_resume(Q,c,op,t,rho,inner_tol,maxiter,verbose,info,state)
        solves++; total=total+info[2]; warm=warm+state.reused; stats=stats+info[12..15]
        T=colshape(t,K)
        slopes=J(n,K,0); slopes[,free]=T[,2..K]
        slopes[,pivot]=(J(n,1,1)-slopes[,free]*a[free]')/a[pivot]
        monos=max((0\-vec(slopes)))
        slopes=slopes:/scales
        coeff=(T[,1],slopes)
        mono=max((0\-vec(slopes)))
        normerr=max(abs(slopes*g':-1))
        d=rowsum(W:*coeff,1)
        rss=sum(d:^2,1)
        newpairs=_cnlsddf_scan(W,coeff,added,tol,viol)
        hist=hist\(cp,rows(pairs),info[1],viol,rows(newpairs),info[2],info[11])
        if (hasmissing((t',d',rss,mono,monos,normerr,viol,info[3..10])) | info[1]<0) {
            reason="numerical_failure"; break
        }
        if (info[1]==1 & info[11]<=tol*.05 & viol<=tol & monos<=tol & normerr<=ntol) {
            conv=1; reason="certified"; break
        }
        if (verbose) printf("{txt}CP %g: active=%g max Afriat=%9.2e normalization=%9.2e\n",cp,rows(pairs),viol,normerr)
        if (cp==maxcut) {
            if (info[1]==0) reason="inner_maxiter"
            else if (normerr>ntol) reason="normalization"
            break
        }
        if (rows(newpairs)) {
            pairs=pairs\newpairs
            for (j=1;j<=rows(newpairs);j++) added[newpairs[j,1],newpairs[j,2]]=1
            nextop=_cnlsddf_pair_operator(B[,2..K],pairs,0,0,1,1,(0,a[free]))
            if (coldstart) state.valid=0
            else expanded=expanded+_cnlsddf_expand_state(op,nextop,t,state)
            op=nextop
            // Intermediate QPs only: never certify at this loose tolerance.
            if (adaptive) inner_tol=max((tol*.05,min((1e-3,.01*viol))))
        }
        else if (inner_tol>tol*.05) inner_tol=tol*.05
        else if (info[1]==1) inner_tol=inner_tol*.1
    }
    obs=selectindex(st_data(.,use))
    st_store(obs,dvar,d); st_store(obs,rvar,-d)
    if (coefvars!="") st_store(obs,tokens(coefvars),coeff)
    summary=(conv,cp,rows(pairs),viol,mono,monos,normerr,rss,total,solves,warm,
        stats,rho,pivot,P,_cnlsddf_op_bytes(op),8*op.ncons*P,mean(d),expanded)
    st_matrix("r(ddf_summary)",summary)
    st_matrix("r(ddf_inner)",info)
    st_matrix("r(ddf_history)",hist)
    st_matrix("r(ddf_gx)",g[1..m]); st_matrix("r(ddf_gy)",g[m+1..K])
    st_local("ddf_reason",reason)
}
end
