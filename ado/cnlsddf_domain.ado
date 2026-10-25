*! version 0.1.0-domain2 02oct2026
program define cnlsddf_domain, rclass sortpreserve
    version 16.0
    syntax , PYthon(string) [GENerate(name) REPLACE]
    if "`e(cmd)'"!="cnlsddf" | "`e(rts)'"!="crs" | "`e(estimator)'"!="cnlsddf_sparse" {
        di as error "cnlsddf_domain requires current sparse CRS CNLS-DDF results"
        exit 301
    }
    if e(converged)!=1 {
        di as error "stage one did not certify convergence"
        exit 498
    }
    if "`e(backend_version)'"!="0.1.0-crs2" | "`e(bridge_version)'"!="0.1.0-dir1" | "`e(source_schema)'"!="ddf_crs_rowhash_v1" {
        di as error "unsupported first-stage version or provenance schema"
        exit 459
    }
    if "`e(coefficient_prefix)'"=="" {
        di as error "refit with coefficients() to retain the fitted planes"
        exit 459
    }
    if !inlist("`c(os)'","MacOSX","Unix") exit 198
    local inputs `e(inputs)'
    local outputs `e(outputs)'
    local distance `e(distance_var)'
    local residual `e(resid_var)'
    local coef `e(coefficient_prefix)'
    local m : word count `inputs'
    local s : word count `outputs'
    if `m'<1 | `s'<1 | e(m)!=`m' | e(s)!=`s' exit 459
    local protected `inputs' `outputs' `distance' `residual' `coef'_alpha
    forvalues j=1/`m' {
        local protected `protected' `coef'_beta`j'
    }
    forvalues j=1/`s' {
        local protected `protected' `coef'_gamma`j'
    }
    if "`e(source_variables)'"!="`protected'" | "`e(source_signature)'"=="" exit 459
    foreach v of local protected {
        capture confirm numeric variable `v', exact
        if _rc exit 459
    }
    tempvar touse
    quietly gen byte `touse'=e(sample)
    quietly count if `touse'
    local n=r(N)
    if `n'!=e(N) | `n'<1 exit 459
    tempname xdirection ydirection
    matrix `xdirection'=e(gx)
    matrix `ydirection'=e(gy)
    if rowsof(`xdirection')!=1 | colsof(`xdirection')!=`m' | rowsof(`ydirection')!=1 | colsof(`ydirection')!=`s' exit 459
    quietly _cnlsddf_signature `protected', sample(`touse') gx(`xdirection') gy(`ydirection')
    if "`r(signature)'"!="`e(source_signature)'" {
        di as error "stage-one sample, values, planes or direction changed; refit before diagnosis"
        exit 459
    }
    if "`generate'"=="" local generate _ddfdom
    local fields rawdist envelope zero lower upper lbinf ubinf step feasible improve observed status endpoint negany negmat negx negy minraw marglo marghi minqty minplane
    foreach kind in raw target {
        forvalues j=1/`m' {
            local fields `fields' `kind'_x`j'
        }
        forvalues j=1/`s' {
            local fields `fields' `kind'_y`j'
        }
    }
    local destinations
    foreach field of local fields {
        local dest `generate'_`field'
        if strlen("`dest'")>32 exit 198
        confirm name `dest'
        if "`replace'"=="" confirm new variable `dest'
        local destinations `destinations' `dest'
    }
    local overlap : list destinations & protected
    if "`overlap'"!="" exit 198
    confirm file `"`python'"'
    quietly findfile cnlsddf_domain.py
    local backend `"`r(fn)'"'
    tempfile input bridge backendlog
    foreach path in python backend input bridge backendlog {
        mata: st_local("unsafe",strofreal(any(strpos(st_local("`path'"),(char(34),char(36),char(96),char(92),char(10),char(13))):>0)))
        if `unsafe' exit 198
    }
    mata: _cnlsddfd_export(st_local("input"),"`protected'","`touse'")
    mata: st_local("gx",invtokens(strtrim(strofreal(st_matrix("`xdirection'"),"%24.17g")))); st_local("gy",invtokens(strtrim(strofreal(st_matrix("`ydirection'"),"%24.17g"))))
    mata: st_local("rss",strtrim(strofreal(st_numscalar("e(rss)"),"%24.17g"))); st_local("mean",strtrim(strofreal(st_numscalar("e(mean_distance)"),"%24.17g")))
    quietly shell "`python'" "`backend'" --input "`input'" --output "`bridge'" --m `m' --s `s' --gx `gx' --gy `gy' --rss `rss' --mean=`mean' >"`backendlog'" 2>&1
    local backend_rc=0
    foreach suffix in csv ok {
        capture confirm file "`bridge'.`suffix'"
        if _rc local backend_rc=498
    }
    if `backend_rc' {
        capture type "`backendlog'"
        capture erase "`bridge'.csv"
        capture erase "`bridge'.ok"
        di as error "domain backend did not complete; no outputs changed"
        exit `backend_rc'
    }
    local nout : word count `destinations'
    local temporary
    forvalues j=1/`nout' {
        tempvar out`j'
        quietly gen double `out`j''=.
        local temporary `temporary' `out`j''
    }
    tempname summary
    capture noisily mata: _cnlsddfd_import(st_local("bridge"),"`fields'","`temporary'","`touse'",`n',`m'+`s',"`summary'")
    local backend_rc=_rc
    capture erase "`bridge'.csv"
    capture erase "`bridge'.ok"
    if `backend_rc' exit `backend_rc'
    // A later output failure restores data, labels, types and characteristics.
    preserve
    _cnlsddf_domain_commit, destinations(`destinations') temporary(`temporary') `replace'
    label variable `generate'_status "Domain: 0 empty, 1 zero surface, 2 input cap, 3 near boundary"
    label variable `generate'_endpoint "Endpoint: 0 none, 1 zero surface, 2 input cap, 3 both"
    label variable `generate'_step "Signed physical ray endpoint; missing if ray empty"
    label variable `generate'_feasible "Nonnegative-quantity ray intersection exists"
    label variable `generate'_improve "Ray permits a nonnegative step in stored technology"
    label variable `generate'_lbinf "Lower bound is negative infinity"
    label variable `generate'_ubinf "Upper bound is positive infinity"
    forvalues j=1/`nout' {
        local dest : word `j' of `destinations'
        char `dest'[cnlsddf_domain] "0.1.0-domain2; `e(source_signature)'"
    }
    local names N empty near_boundary input_cap raw_negative negative_step
    local j=0
    foreach name of local names {
        local ++j
        return scalar `name'=`summary'[1,`j']
    }
    return local cmd cnlsddf_domain
    return local version "0.1.0-domain2"
    return local source_status "ddf_crs_rows_verified"
    return local generate "`generate'"
    di as text "CRS domain diagnostics: N=" as result `summary'[1,1] as text "  empty rays=" as result `summary'[1,2] as text "  near boundary=" as result `summary'[1,3]
    restore, not
end

program define _cnlsddf_domain_commit
    version 16.0
    syntax , DESTinations(namelist) TEMPOrary(namelist) [REPLACE]
    local n : word count `destinations'
    forvalues j=1/`n' {
        local dest : word `j' of `destinations'
        local tmp : word `j' of `temporary'
        if "`replace'"!="" {
            capture confirm variable `dest', exact
            if !_rc drop `dest'
        }
        rename `tmp' `dest'
    }
end

version 16.0
mata:
mata set matastrict on
void _cnlsddfd_export(string scalar path,string scalar vars,string scalar sample)
{
    real matrix z
    real scalar f,i
    z=st_data(.,tokens(vars),sample)
    if (hasmissing(z)) _error(459,"Missing stage-one values")
    f=fopen(path,"w"); fput(f,invtokens(tokens(vars),","))
    for(i=1;i<=rows(z);i++) fput(f,invtokens(strtrim(strofreal(z[i,.],"%24.17g")),","))
    fclose(f)
}
void _cnlsddfd_import(string scalar stem,string scalar fields,string scalar destinations,string scalar sample,real scalar n,real scalar k,string scalar summary)
{
    real scalar f,i,j
    real matrix z,expected
    real colvector feasible
    string matrix line
    f=fopen(stem+".ok","r"); line=fget(f); fclose(f)
    if (line!="cnlsddf_domain_v1") _error(498,"Invalid completion marker")
    f=fopen(stem+".csv","r"); line=fget(f)
    if (line!=invtokens(tokens(fields),",")) _error(498,"Invalid domain header")
    z=J(0,cols(tokens(fields)),.)
    while((line=fget(f))!=J(0,0,"")) {
        if (cols(tokens(subinstr(line,","," ")))!=cols(z)) _error(498,"Invalid domain columns")
        z=z\strtoreal(tokens(subinstr(line,","," ")))
    }
    fclose(f)
    if (rows(z)!=n | cols(z)!=22+2*k) _error(498,"Invalid domain dimensions")
    // Exactly specified missing values: unbounded limits/margins and absent endpoints.
    expected=J(n,cols(z),0)
    expected[,4]=z[,6]; expected[,5]=z[,7]
    expected[,19]=z[,6]; expected[,20]=z[,7]
    feasible=z[,9]
    expected[,8]=expected[,21]=expected[,22]=1:-feasible
    for(j=23+k;j<=22+2*k;j++) expected[,j]=1:-feasible
    if (any((z:>=.):!=expected)) _error(498,"Invalid domain missingness")
    for(j=6;j<=17;j++) {
        if (any(j:==(8,12,13))) continue
        if (any((z[,j]:!=0):&(z[,j]:!=1))) _error(498,"Invalid domain flag")
    }
    if (any(z[,12]:<0) | any(z[,12]:>3) | any(z[,12]:!=floor(z[,12])) | any(z[,13]:<0) | any(z[,13]:>3) | any(z[,13]:!=floor(z[,13]))) _error(498,"Invalid domain state")
    if (any((z[,13]:==0):!=(!feasible)) | any((z[,12]:==0):&(feasible:==1)) | any((z[,10]:==1):&((feasible:==0):|(z[,8]:<0)))) _error(498,"Inconsistent domain state")
    st_store(.,tokens(destinations),sample,z)
    st_matrix(summary,(n,sum(!feasible),sum(z[,12]:==3),sum(z[,12]:==2),sum(z[,15]),sum(z[,8]:<0)))
}
end
