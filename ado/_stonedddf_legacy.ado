*! version 0.1.0-stat1 01oct2026
*! Distance-aware normal/half-normal MM, PSL and QLE postestimation for certified cnlsddf.
program define _stonedddf_legacy, rclass sortpreserve
    version 16.0
    syntax , [GENerate(name) TARgets(name) NOISEAdjusted(name) REPLACE ///
        METHOD(string) MAXiter(integer 200) LTOLerance(real 1e-8)]
    if "`method'"=="" local method mm
    local method=lower("`method'")
    if !inlist("`method'","mm","psl","qle") {
        di as error "method() must be mm, psl or qle"
        exit 198
    }
    if missing(`maxiter') | `maxiter'<1 | missing(`ltolerance') | `ltolerance'<=0 | `ltolerance'>1e-3 {
        di as error "maxiter() must be positive; ltolerance() must be in (0, .001]"
        exit 198
    }
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
    tempvar touse
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
    local destinations `generate'_u `generate'_v `generate'_distance
    foreach prefix in targets noiseadjusted {
        if "``prefix''"!="" {
            forvalues j=1/`m' {
                local destinations `destinations' ``prefix''_x`j'
            }
            forvalues j=1/`s' {
                local destinations `destinations' ``prefix''_y`j'
            }
            local destinations `destinations' ``prefix''_negative
        }
    }
    local unique : list uniq destinations
    local nout : word count `destinations'
    local nuniq : word count `unique'
    local overlap : list destinations & protected
    if `nout'!=`nuniq' | "`overlap'"!="" {
        di as error "output names must be distinct and cannot replace stage-one variables"
        exit 198
    }
    foreach v of local destinations {
        confirm name `v'
        if strlen("`v'")>32 error 198
        if "`replace'"=="" confirm new variable `v'
    }
    local temporary
    forvalues j=1/`nout' {
        tempvar out`j'
        quietly gen double `out`j''=.
        local temporary `temporary' `out`j''
    }
    tempname results likelihood profile
    mata: _stonedddf_run("`touse'", "`inputs'", "`outputs'", "`distance'", "`residual'", ///
        "`temporary'", ("`targets'"!=""), ("`noiseadjusted'"!=""), "`results'", "`method'", `ltolerance', `maxiter', "`likelihood'", "`profile'")
    // Automatic restore on any later error also restores types, labels and characteristics.
    preserve
    // No destination is changed until all parameters and generated values are finite.
    forvalues j=1/`nout' {
        local dest : word `j' of `destinations'
        local tmp : word `j' of `temporary'
        if "`replace'"!="" {
            capture confirm variable `dest', exact
            if !_rc drop `dest'
        }
        rename `tmp' `dest'
    }
    label variable `generate'_u "Directional E(u | v-u), `method'"
    label variable `generate'_v "Directional E(v | v-u), `method'"
    label variable `generate'_distance "Signed distance to shifted StoNED frontier, `method'"
    foreach prefix in targets noiseadjusted {
        if "``prefix''"!="" {
            label variable ``prefix''_negative "Any negative projected quantity (no clipping)"
        }
    }
    local names N sigma_u sigma_v mu mean_distance frontier_shift m2 m3 boundary target_negative noise_negative
    local j=0
    foreach name of local names {
        local ++j
        return scalar `name'=`results'[1,`j']
    }
    return local source_status "ddf_rows_verified"
    return local source_schema "`e(source_schema)'"
    return local cmd stonedddf
    if "`method'"!="mm" {
        local lnames ll converged search_iterations evaluations noise_ratio ll_gain boundary_tolerance normal_ll candidate_ll
        local j=0
        foreach name of local lnames {
            local ++j
            return scalar `name'=`likelihood'[1,`j']
        }
        matrix colnames `profile'=tau sigma_u sigma_v ll
        return matrix likelihood_profile=`profile'
        return scalar ltolerance=`ltolerance'
        return scalar noise_floor=1e-5
        return local selection "`sd_status'"
    }
    return local method "`method'"
    return local distribution normal_halfnormal
    return local status "`sd_status'"
    return local generate "`generate'"
    return local targets "`targets'"
    return local noiseadjusted "`noiseadjusted'"
    di as text "Directional StoNED (`method'): sigma_u=" as result %10.6g `results'[1,2] ///
        as text "  sigma_v=" as result %10.6g `results'[1,3] as text "  (`sd_status')"
    if `results'[1,9]>0 di as text "Boundary fit: this is not evidence that all units are efficient."
    if `results'[1,10]>0 & `results'[1,10]<. di as text "Negative frontier target quantities: " as result `results'[1,10]
    if `results'[1,11]>0 & `results'[1,11]<. di as text "Negative noise-adjusted quantities: " as result `results'[1,11]
    restore, not
end

version 16.0
mata:
mata set matastrict on

// E[Z | Z>=0], Z~N(a,1). The continued fraction evaluates the SMALL
// difference phi(t)/Phi(-t)-t directly, avoiding tail cancellation/underflow.
real scalar _stonedddf_tmean(real scalar a)
{
    real scalar k, tail, t
    if (missing(a)) return(.)
    if (a>=-8) return(a+exp(lnnormalden(a)-lnnormal(a)))
    t=-a
    tail=0
    for (k=100;k>=2;k--) tail=k/(t+tail)
    return(1/(t+tail))
}

real colvector _stonedddf_posterior(real colvector z, real scalar su, real scalar sv)
{
    real scalar scale, sig, ratio, ss, i, a
    real colvector u
    if (su<0 | sv<=0 | hasmissing((su,sv,z'))) _error(498)
    if (su==0) return(J(rows(z),1,0))
    scale=max((su,sv))
    sig=scale*sqrt((su/scale)^2+(sv/scale)^2)
    ratio=su/sig
    // Compute su*sv/sig without underflowing either scale ratio first.
    ss=min((su,sv))*(scale/sig)
    u=J(rows(z),1,.)
    for (i=1;i<=rows(z);i++) {
        a=-(z[i]/sv)*ratio
        if (a>8 & a<.) u[i]=-z[i]*ratio^2 + ss*exp(lnnormalden(a)-lnnormal(a))
        else u[i]=ss*_stonedddf_tmean(a)
    }
    if (hasmissing(u) | min(u)<0) _error(498)
    return(u)
}

// Return original-unit moments, scales and explicit identification status.
// Population central moments (both denominators N), with scaled arithmetic.
real rowvector _stonedddf_mm(real colvector d)
{
    real scalar scale, dm, v2, v3, su, sv2, sv, mu, boundary, m2, m3
    real colvector r
    if (rows(d)<3 | hasmissing(d)) _error(498)
    scale=max(abs(d))
    if (scale==0) {
        errprintf("Zero residual variance: no identified normal/half-normal decomposition.\n")
        _error(498)
    }
    dm=mean(d/scale)
    r=dm:-d/scale
    v2=mean(r:^2)
    v3=mean(r:^3)
    if (v2==0) {
        errprintf("Zero residual variance: no identified normal/half-normal decomposition.\n")
        _error(498)
    }
    boundary=0
    if (v3>=0) {
        su=0
        boundary=(v3>0 ? 1 : 2)
    }
    else su=(-v3/((4/pi()-1)*sqrt(2/pi())))^(1/3)
    sv2=v2-su^2*(1-2/pi())
    if (sv2<=0) {
        errprintf("MM implies nonpositive noise variance; no interior normal/half-normal decomposition.\n")
        errprintf("No outputs changed. Reconsider the model or try method(psl) or method(qle).\n")
        _error(498)
    }
    sv=sqrt(sv2)*scale
    su=su*scale
    dm=dm*scale
    mu=su*sqrt(2/pi())
    m2=(v2*scale)*scale
    m3=((v3*scale)*scale)*scale
    if (hasmissing((su,sv,mu,dm,m2,m3))) _error(498)
    return((su,sv,mu,dm,mu-dm,m2,m3,boundary))
}

// Centered normal/half-normal likelihood. Stable log Phi follows the
// reviewed GSFA approach; no probability floor and no external Mata globals.
real colvector _stonedddf_mills(real colvector a)
{
    real colvector out, ids, t, tail
    real scalar k
    out=exp(lnnormalden(a)-lnnormal(a))
    ids=selectindex(a:<-8)
    if (rows(ids)) {
        t=-a[ids]
        tail=J(rows(t),1,0)
        for (k=60;k>=2;k--) tail=k:/(t+tail)
        out[ids]=t+1:/(t+tail)
    }
    return(out)
}

// r has mean zero and variance one. For fixed tau=sin(angle), PSL is
// strictly convex in inverse total scale k. Its root lies in [cos(angle),1].
// Return mean negative log likelihood, k, inner iterations.
real rowvector _stonedddf_profile(real scalar angle, real colvector r,
    real scalar psl, real scalar tol, real scalar maxiter, real scalar evaluations)
{
    real scalar tau, lam, k, lo, hi, j, score, f
    real colvector w, a
    tau=sin(angle)
    lam=tan(angle)
    j=0
    if (!psl) k=sqrt(1-2/pi()*tau^2)
    else if (angle==0) k=1
    else {
        lo=cos(angle)
        hi=1
        for (j=1;j<=maxiter;j++) {
            k=(lo+hi)/2
            w=r*k:-sqrt(2/pi())*tau
            a=-lam*w
            score=-1/k+mean(r:*(w+lam*_stonedddf_mills(a)))
            evaluations++
            if (missing(score)) _error(498)
            if (score>0) hi=k
            else lo=k
            if ((hi-lo)<=tol*k) break
        }
        if (j>maxiter) {
            errprintf("PSL scale profile exhausted maxiter(); no outputs changed.\n")
            _error(430)
        }
        k=(lo+hi)/2
    }
    w=r*k:-sqrt(2/pi())*tau
    f=-log(k)-log(2)+.5*log(2*pi())+mean(.5*w:^2-lnnormal(-lam*w))
    evaluations++
    if (missing(f)) _error(498)
    return((f,k,j))
}

// Golden-section refinement of each grid-identified basin, with a width
// certificate. Endpoint candidates are evaluated independently, never inferred
// from an optimizer return code or a nearly flat likelihood value.
real rowvector _stonedddf_refine(real scalar left, real scalar right,
    real colvector r, real scalar psl, real scalar tol, real scalar maxiter,
    real scalar evaluations)
{
    real scalar ratio, x1, x2, j
    real rowvector p1, p2, final
    ratio=(sqrt(5)-1)/2
    x1=right-ratio*(right-left)
    x2=left+ratio*(right-left)
    p1=_stonedddf_profile(x1,r,psl,tol,maxiter,evaluations)
    p2=_stonedddf_profile(x2,r,psl,tol,maxiter,evaluations)
    for (j=1;j<=maxiter;j++) {
        if (right-left<=tol) break
        if (p1[1]<p2[1]) {
            right=x2
            x2=x1
            p2=p1
            x1=right-ratio*(right-left)
            p1=_stonedddf_profile(x1,r,psl,tol,maxiter,evaluations)
        }
        else {
            left=x1
            x1=x2
            p1=p2
            x2=left+ratio*(right-left)
            p2=_stonedddf_profile(x2,r,psl,tol,maxiter,evaluations)
        }
    }
    if (right-left>tol) {
        errprintf("Likelihood angle search exhausted maxiter(); no outputs changed.\n")
        _error(430)
    }
    x1=(right+left)/2
    final=_stonedddf_profile(x1,r,psl,tol,maxiter,evaluations)
    return((x1,final[1..2],min((j-1,maxiter))))
}

real rowvector _stonedddf_likelihood(real colvector d, string scalar method,
    real scalar tol, real scalar maxiter, real rowvector diag, real matrix grid)
{
    real colvector r, agrid
    real scalar ds, dm, sd, m2, m3, top, i, evaluations, iters, bestf, bestangle, bestk
    real scalar su, sv, mu, boundary, k, f, psl, gain, candidatef, boundarytol
    real rowvector p, refined
    real matrix raw
    if (rows(d)<3 | hasmissing(d)) _error(498)
    ds=max(abs(d))
    if (ds==0) {
        errprintf("Zero residual variance: likelihood scale has no finite maximizer.\n")
        _error(498)
    }
    dm=mean(d/ds)*ds
    r=dm/ds:-d/ds
    sd=sqrt(mean(r:^2))*ds
    if (sd<=0 | missing(sd)) {
        errprintf("Zero/invalid residual variance: likelihood is not identified.\n")
        _error(498)
    }
    r=r*(ds/sd)
    m2=sd^2
    m3=((mean(r:^3)*sd)*sd)*sd
    if (hasmissing((m2,m3))) _error(498)
    psl=(method=="psl")
    top=acos(1e-5)
    agrid=(0::64)*(top/64)
    raw=J(rows(agrid),2,.)
    evaluations=iters=0
    for (i=1;i<=rows(agrid);i++) {
        p=_stonedddf_profile(agrid[i],r,psl,tol,maxiter,evaluations)
        raw[i,.]=p[1..2]
    }
    bestangle=0
    bestf=raw[1,1]
    bestk=1
    for (i=1;i<=rows(agrid);i++) {
        if (raw[i,1]<bestf) {
            bestf=raw[i,1]
            bestangle=agrid[i]
            bestk=raw[i,2]
        }
        if (i==1 | i==rows(agrid)) continue
        if (i==2 | i==rows(agrid)-1 | (raw[i,1]<=raw[i-1,1] & raw[i,1]<=raw[i+1,1])) {
            refined=_stonedddf_refine(agrid[i-1],agrid[i+1],r,psl,tol,maxiter,evaluations)
            iters=iters+refined[4]
            if (refined[2]<bestf) {
                bestangle=refined[1]
                bestf=refined[2]
                bestk=refined[3]
            }
        }
    }
    // Statistical selection uses a fixed per-observation numerical tie rule.
    // ltolerance() controls search widths, not the definition of a normal boundary.
    boundarytol=1e-8
    candidatef=bestf
    gain=raw[1,1]-candidatef
    boundary=0
    if (gain<=boundarytol) {
        bestangle=0
        bestf=raw[1,1]
        bestk=1
        boundary=4
    }
    else if (raw[rows(raw),1]-bestf<=boundarytol | top-bestangle<=4*tol) {
        errprintf("Likelihood reaches the zero-noise search limit (sigma_v/sigma=1e-5).\n")
        errprintf("The zero-noise support-boundary model is not certified; no outputs changed.\n")
        _error(498)
    }
    su=sd*sin(bestangle)/bestk
    sv=sd*cos(bestangle)/bestk
    mu=su*sqrt(2/pi())
    diag=(-rows(d)*(bestf+log(sd)),1,iters,evaluations,cos(bestangle),
        gain,boundarytol,-rows(d)*(raw[1,1]+log(sd)),-rows(d)*(candidatef+log(sd)))
    grid=J(rows(agrid),4,.)
    for (i=1;i<=rows(agrid);i++) {
        k=raw[i,2]
        f=raw[i,1]
        grid[i,.]=(sin(agrid[i]),sd*sin(agrid[i])/k,sd*cos(agrid[i])/k,-rows(d)*(f+log(sd)))
    }
    if (hasmissing((su,sv,mu,dm,m2,m3,diag))) _error(498)
    return((su,sv,mu,dm,mu-dm,m2,m3,boundary))
}

void _stonedddf_run(string scalar sample, string scalar inputs, string scalar outputs,
    string scalar distance, string scalar residual, string scalar destinations,
    real scalar targets, real scalar noiseadjusted, string scalar result,
    string scalar method, real scalar tol, real scalar maxiter, string scalar diagnostics, string scalar gridname)
{
    real colvector d, r, D, z, u, v, bad
    real matrix X, Y, out, projected, grid
    real rowvector gx, gy, mm, diag
    real scalar n, rss, tn, nn, ds, gap, meanold
    string rowvector statuses
    d=st_data(.,distance,sample)
    r=st_data(.,residual,sample)
    n=rows(d)
    ds=max((1,max(abs(d))))
    gap=max(abs(d/ds+r/ds))
    rss=sum(d:^2,1)
    meanold=st_numscalar("e(mean_distance)")
    if (hasmissing((gap,rss,meanold)) | gap>1e-12 |
        abs(rss-st_numscalar("e(rss)"))>1e-12*max((1,rss)) |
        abs(mean(d/ds)-meanold/ds)>1e-12) {
        errprintf("Stage-one distance/residual diagnostics changed; refit cnlsddf.\n")
        _error(459)
    }
    if (method=="mm") mm=_stonedddf_mm(d)
    else {
        mm=_stonedddf_likelihood(d,method,tol,maxiter,diag,grid)
        st_matrix(diagnostics,diag)
        st_matrix(gridname,grid)
    }
    D=(d:-mm[4]):+mm[3]
    z=-D
    u=_stonedddf_posterior(z,mm[1],mm[2])
    v=z+u
    out=(u,v,D)
    tn=nn=.
    X=st_data(.,tokens(inputs),sample)
    Y=st_data(.,tokens(outputs),sample)
    gx=st_matrix("e(gx)")
    gy=st_matrix("e(gy)")
    if (targets) {
        projected=(X-D*gx,Y+D*gy)
        if (hasmissing(projected)) _error(498)
        bad=rowmin(projected):<0
        tn=sum(bad)
        out=(out,projected,bad)
    }
    if (noiseadjusted) {
        projected=(X+v*gx,Y-v*gy)
        if (hasmissing(projected)) _error(498)
        bad=rowmin(projected):<0
        nn=sum(bad)
        out=(out,projected,bad)
    }
    if (hasmissing(out) | hasmissing(mm)) _error(498)
    st_store(.,tokens(destinations),sample,out)
    st_matrix(result,(n,mm,tn,nn))
    statuses=("interior","wrong_skew","zero_skew","zero_variance","normal_boundary")
    st_local("sd_status",statuses[mm[8]+1])
}
end
