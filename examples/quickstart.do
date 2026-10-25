version 16
// In a fresh Stata session:
// do "/absolute/path/CNLS-StoNED-Stata/examples/quickstart.do" "/absolute/path/CNLS-StoNED-Stata"
args root
if `"`root'"'=="" {
    display as error "Supply the repository's absolute path as the first argument."
    exit 198
}
adopath ++ `"`root'/ado"'
preserve
clear
// Six artificial observations: a syntax/numerical example, not an inference study.
input double(x y)
1 6
1 5.8
1 5.2
1 4.9
1 4.5
1 3.6
end
convnls y, inputs(x) frontier(production) generate(ehat) fitted(fhat) ///
    init(ele) tolerance(1e-10) maxiter(100000) maxcut(300)
assert e(converged)==1
stoned ehat, outcome(y) generate(reference_score) inefficiency(reference_u)
assert e(parameter_valid)==1 & e(score_valid)==1 & e(estimate_valid)==1
stoned ehat, implementation(legacy) generate(bc_score) inefficiency(legacy_u)
// These scores have different definitions; their equality is not expected.
list y ehat reference_score bc_score, noobs
restore
display "QUICKSTART PASSED"
