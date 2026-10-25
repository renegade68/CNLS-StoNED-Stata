*! version 0.1.0-safety1 01oct2026
*! Order-invariant full-row checksums; ordinary change detection, not authentication.
program define _cnlsddf_signature, rclass
    version 16.0
    syntax varlist(numeric), SAMPLE(varname numeric) GX(name) GY(name)
    tempvar forward reverse
    quietly gen double `forward'=.
    quietly gen double `reverse'=.
    mata: _cnlsddf_rowhash("`varlist'", "`sample'", "`gx'", "`gy'", "`forward' `reverse'")
    quietly _datasignature `forward' `reverse' if `sample', nonames
    return local signature "`r(datasignature)'"
end
version 16.0
mata:
mata set matastrict on
void _cnlsddf_rowhash(string scalar vars, string scalar sample,
    string scalar gx, string scalar gy, string scalar destination)
{
    real matrix z, out
    real rowvector row, direction
    real scalar i
    z=st_data(.,tokens(vars),sample)
    direction=(st_matrix(gx),st_matrix(gy))
    out=J(rows(z),2,.)
    for (i=1;i<=rows(z);i++) {
        row=(z[i,.],direction)
        // Explicit byte order and double values make storage format irrelevant.
        out[i,1]=hash1(row,.,2)
        out[i,2]=hash1(row[cols(row)..1],.,2)
    }
    st_store(.,tokens(destination),sample,out)
}
end
