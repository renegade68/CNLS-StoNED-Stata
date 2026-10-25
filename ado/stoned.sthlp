{smcl}
{* Public packaging documentation revision 25oct2026; estimator unchanged.}{...}
{* *! version 0.4.3-integrated1-doc1 24oct2026}{...}
{title:stoned -- pyStoNED reproduction, raw diagnostics, and legacy estimates}

{pstd}
The reference is the unmodified JSS pyStoNED 0.7.3 release.
Reference compatibility does not establish statistical validity. Requires Stata 16; the
default MOM and legacy MM/PSL do not call Python or an external solver.
QLE uses a source-checked external Python process with the published optimizer.

{title:Syntax}
{p 8 17 2}
{cmd:stoned} {it:residual} {ifin} [{cmd:,}
{opt implementation(pystoned|legacy)} {opt method(mm|qle|psl)} {opt invalid(refuse|raw)}
{opt frontier(production|cost)} {opt cet(addi|mult)} {opt outcome(varname)}
{opt pythonexec(path)} {opt diagnose}
{opt generate(newvar)} {opt inefficiency(newvar)} {opt verbose}]

{title:Default and reference scope}
{pstd}
Defaults are {cmd:implementation(pystoned) method(mm) frontier(production)
cet(addi)}. ADDI requires {cmd:outcome(y)}: the original outcome on the
same estimation sample, in the same units as the residual. Missing or zero y
is refused without silently reducing the residual sample.
{cmd:method(mom)} is an alias for {cmd:method(mm)}.
The default {cmd:invalid(refuse)} preserves the prior safety behavior. With
{cmd:implementation(pystoned) invalid(raw)}, the command retains raw public
parameters and optimizer state in e(), returns successfully with
{cmd:e(parameter_valid)=0} for inadmissible parameters, and marks scores missing
when the reference calculation is nonfinite. This mode is for reproduction
diagnostics; its scores are not certified estimates. It does not silently
repair parameters or choose another estimator.

{pstd}
{cmd:e(parameter_valid)} tests finite nonnegative sigma_u and finite positive
sigma_v. {cmd:e(score_valid)} tests finite scores/conditional quantities and
nonnegative conditional U. {cmd:e(estimate_valid)} combines these checks with
QLE optimizer success. Neither flag is evidence of statistical accuracy.
In particular, a nonconverged QLE can have parameter_valid=1 and score_valid=1
while estimate_valid=0. Raw variables are labelled as reference quantities.
{cmd:e(score_status)=reference_nonfinite_score} identifies nonfinite score
cells separately from {cmd:reference_negative_u}.

{pstd}
The compatibility MOM uses centered empirical moments with divisor N.
It follows the fixed reference's wrong-skew convention: production M3 is set
to zero; cost M3 is replaced with 0.00001. This replacement is not unit
invariant. Its use is visible in {cmd:e(m3_raw)} and {cmd:e(m3_used)}.
QLE calls the unchanged public implementation with one unconstrained lambda,
raw mean squared residuals and initial lambda=1. Specify {cmd:pythonexec(path)}
to a Python executable containing pyStoNED 0.7.3, NumPy 2.0.2 and SciPy 1.13.1.
The default executable name is python3; it must satisfy these checks. All 39
public Python source hashes are checked before and after each call. No fallback
is used. The external bridge was tested on macOS; the external launch path has
not been validated on Windows. No Stata Python shared-library configuration
is required and the command does not change persistent Python settings.

{pstd}
Let U be the reference conditional-inefficiency calculation. The reference
score includes the published 1e-6 denominator adjustment. ADDI returns
(y-U)/y for production and (y+U)/y for cost. MULT returns exp(-U) and
exp(U), respectively. These are not BC E[exp(-u)|z] scores and are not
clipped to [0,1]. Negative conditional U is preserved for reproducibility
and flagged with {cmd:e(score_valid)=0}; such output is not certified as
valid conditional inefficiency. Negative/nonfinite sigma_u, nonpositive
sigma_v or nonconverged QLE is refused under the default policy. Under
{cmd:invalid(raw)}, raw parameters/status are returned separately, with
{cmd:e(parameter_valid)}, {cmd:e(score_valid)} and {cmd:e(invalid_reason)}
distinguishing parameter and score validity. NaN/infinite Python values are
represented by Stata missing scalars plus literal text in
{cmd:e(raw_sigma_u_text)}, {cmd:e(raw_sigma_v_text)} and
{cmd:e(raw_lambda_text)}; they are not claimed as exact numeric Stata values.
QLE also records {cmd:e(qle_optimizer_success)},
{cmd:e(qle_optimizer_status)}, {cmd:e(qle_optimizer_nfev)} and
{cmd:e(qle_optimizer_message)}.
The checked bridge response must contain a complete versioned record. Missing,
truncated, malformed or unavailable bridge output is refused even in raw mode.
Only admissible public results count as score agreement.

{pstd}
{cmd:cet(mult)} selects the public score transformation only. It does not
refit a multiplicative first-stage CNLS model. Input residuals must come
from the model intended by the user. KDE, DDF-specific projected outcomes,
and the public frontier-envelope method are outside this command's new scope.

{title:Legacy implementation}
{pstd}
{cmd:implementation(legacy)} preserves RC3 sample moments, MM/PSL, stable
conditional formulas and BC scoring. Its method is {cmd:mm} (default) or
{cmd:psl}. Without {cmd:diagnose}, {cmd:cet()} and {cmd:outcome()}
are rejected. With {cmd:diagnose}, these options specify only the public
comparison mapping: diagnostic {cmd:cet(addi)} requires {cmd:outcome(y)},
whereas explicit {cmd:cet(mult)} omits {cmd:outcome()}.
The selected legacy fit and BC scores retain their original definitions.
{cmd:invalid(raw)} remains unavailable with {cmd:implementation(legacy)}.
This is the independent prior Stata method, not a corrected pyStoNED mode.
Use this option explicitly when reproducing earlier manuscript tables and
figures. {help stoned_legacy:Legacy formula documentation} describes that
implementation; its old default wording applies only within legacy mode.

{title:Examples}
{phang2}{cmd:convnls y, inputs(x1 x2) generate(ehat)}{p_end}
{phang2}{cmd:stoned ehat, outcome(y) generate(te)}{p_end}
{phang2}{cmd:stoned ehat, outcome(y) method(qle) pythonexec("/path/to/python") generate(te_qle)}{p_end}
{phang2}{cmd:stoned ehat, outcome(y) method(qle) invalid(raw) pythonexec("/path/to/python") generate(ref_score)}{p_end}
{phang2}{cmd:display e(parameter_valid) " " e(score_valid) " " e(estimate_valid)}{p_end}
{phang2}{cmd:stoned ehat, implementation(legacy) generate(te_bc)}{p_end}
{phang2}{cmd:stoned ehat, implementation(legacy) method(psl) generate(te_psl)}{p_end}
{phang2}{cmd:stoned ehat, implementation(legacy) outcome(y) diagnose pythonexec("/path/to/python") generate(te_bc_diag)}{p_end}
{phang2}{cmd:stoned ehat, implementation(legacy) cet(mult) diagnose pythonexec("/path/to/python") generate(te_bc_mult)}{p_end}

{title:Stored results and safety}
{pstd}
In addition to sigma_u, sigma_v, lambda, eu, N, m2, m3 and score summaries,
e() identifies implementation, reference_version, score_definition, cet,
outcome, score_valid, m3_raw and m3_used. QLE also stores qle_converged.
Residual provenance, frontier consistency, output-name protection and
rollback follow the RC3 contract. Input and output variables must be distinct.
Legacy-only PSL results retain their earlier meaning. A valid score status is
not a first-stage convergence certificate or a statistical confidence statement.

{title:Comparison diagnostics}
{pstd}
{cmd:diagnose} restores the established four-by-seventeen {cmd:e(diagnose)}
comparison table: public_mm, public_qle, legacy_mm, legacy_psl. Its valid column
is parameter validity only. Score metrics are evaluated only for admissible
parameters; unavailable metrics remain missing even when raw outputs are retained.
{cmd:invalid(raw)} also applies to the two public diagnostic rows, so rc=0
can coexist with valid=0. The new {cmd:e(diagnose_validity)} table has two public
rows and parameter_valid, score_valid, estimate_valid columns. With refuse,
the QLE parameter probe does not certify the requested score mapping; score_valid
is missing and estimate_valid is zero on a refused probe, otherwise missing.
{cmd:e(diagnose_invalid_policy)} records the policy. Diagnostics never select or
repair the requested fit. Final scalar estimates remain in e(); directional
postestimation uses r() and preserves the first-stage e().
{pstd}
The original diagnostic columns are rc, sigma_u, sigma_v, lambda, valid,
converged, nll_mean, moment_gap, U_mean, BC_mean, mapped_mean, literal_mean,
max_U_gap, delta_zero, literal_U_min, negative_U_n, flat_zero. BC, stable mapped
and literal public scores have different definitions; lower NLL is not an
accuracy comparison or an automatic model-selection rule. Boundary flatness
is a numerical screen, not an inferential test. Legacy PSL reports its effective
method and selection separately. Raw residual levels are not recentered.
