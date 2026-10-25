# CNLS and StoNED for Stata

Software accompanying Choonjoo Lee, *convnls and stoned: Stochastic Non-smooth
Envelopment of Data in Stata* (accepted, Stata Journal), with subsequent author
updates. See [changes](CHANGELOG.md) before updating an older installation.

| Commands | Scope |
| --- | --- |
| [`convnls`](ado/convnls.sthlp), [`stoned`](ado/stoned.sthlp) | Scalar CNLS and noise/inefficiency decomposition; production and cost |
| [`cnlsddf`](ado/cnlsddf.sthlp), [`stonedddf`](ado/stonedddf.sthlp) | Experimental multiple-output directional models; production VRS |
| [`cnlsddf_sparse`](ado/cnlsddf_sparse.sthlp), [`cnlsddf_domain`](ado/cnlsddf_domain.sthlp) | Optional Python engine; CRS first-stage and geometry diagnostics |

## Install and use

Requires **Stata 16+**. Download or clone this repository, keep **all files in
`ado/` together**, and start a fresh Stata session after updating:

```stata
adopath ++ "/your/path/CNLS-StoNED-Stata/ado"
which convnls
which stoned
help stoned
```

For a production model with outcome `y` and inputs `x1 x2` already loaded:

```stata
convnls y, inputs(x1 x2) frontier(production) generate(ehat)
assert e(converged)==1
stoned ehat, outcome(y) generate(score)
display e(parameter_valid) " " e(score_valid) " " e(estimate_valid)
```

[`examples/quickstart.do`](examples/quickstart.do) provides a self-contained
synthetic example. Choose unused output names; scalar commands can replace
existing output variables.

`stoned` now defaults to **pyStoNED 0.7.3 compatibility MOM**, with
`outcome(y)` in the same units as the residual. For the earlier Stata estimator
family and BC scores, specify `implementation(legacy)` explicitly;
`method(psl)` is available in that mode. The two modes return different quantities.

Native MOM and legacy methods need no Python. Public-reference QLE (including
its `diagnose` row) and the sparse engine require separate optional environments;
see [dependencies](DEPENDENCIES.txt).

## Interpretation

Inspect first-stage convergence and second-stage validity/boundary diagnostics.
Reference scores are not necessarily between zero and one. Compatibility and
numerical validity flags do not establish statistical accuracy. Context effects
are point estimates; validated standard errors, tests and confidence intervals
are not provided.

DDF remains a research extension. Its default reference path uses raw distance
and coefficient-dependent ratios with unresolved interpretation/identification
limitations. Use explicit `implementation(legacy)` for directional quantities
and targets; see [DDF help](ado/stonedddf.sthlp). CRS StoNED decomposition is
not supported.

## Attribution and feedback

[AFL-3.0](LICENSE); see [software notice](LICENSE_CODE.txt). Cite the article and
method references in the help files. Journal acceptance does not validate later
extensions. Final citation/DOI will be added when available.

For [issues](https://github.com/renegade68/CNLS-StoNED-Stata/issues), include the
Git commit, Stata version, command, diagnostics and a small reproducible example.

Contact: Choonjoo Lee — bloom.rampike@gmail.com
