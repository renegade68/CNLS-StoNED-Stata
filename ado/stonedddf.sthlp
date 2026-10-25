{smcl}
{* Public packaging documentation revision 25oct2026; estimator unchanged.}{...}
{* *! help revision ddf-interpretation1 25oct2026; executable 0.2.1-rawdiag1 unchanged}{...}
{title:stonedddf -- Fixed pyStoNED DDF reference and explicit legacy option}

{p 4 4 2}
Research extension: production VRS, desirable outputs, common direction.
Run after the included {cmd:cnlsddf} or a certified VRS {cmd:cnlsddf_sparse} fit.
Keep all files in ado/ together and start a fresh Stata session after updating.
The full original sample and same-fit coefficients are required for public scoring.
This is not a new first-stage solver or the journal scalar stoned command.

{title:Syntax}
{p 8 12 2}
{cmd:stonedddf}, [{opt implementation(pystoned|legacy)} {opt method(mm|qle|psl)}
{opt invalid(refuse|raw)} {opt generate(prefix)} {opt replace} {opt diagnose} {opt pythonexec(path)}
{opt targets(prefix)} {opt noiseadjusted(prefix)} {opt maxiter(#)} {opt ltolerance(#)}]

{p 4 4 2}
Default {cmd:implementation(pystoned)} and {cmd:method(mm)} reproduce literal
pyStoNED 0.7.3 MOM/ADDI on d (the first-stage distance), using gamma_i'y_i from
the same saved coefficients. Generates {cmd:prefix_score} and {cmd:prefix_raw_u}.
Default prefix is {cmd:_stonedddf}. The reference ratio is not a directional distance,
Farrell score, BC score, or geometric projection. Negative literal conditional
quantities are retained and flagged; there is no clipping.

{p 4 4 2}
{cmd:method(qle)} uses the unchanged, source-checked public Python method:
pyStoNED 0.7.3, NumPy 2.0.2, SciPy 1.13.1. {cmd:pythonexec()} identifies the local
Python executable (default python3). There is no installation, network request,
automatic fallback or conversion to the legacy optimizer. MOM needs no Python.
{cmd:psl}, {cmd:targets()} and {cmd:noiseadjusted()} require the legacy option.
{cmd:maxiter()} and {cmd:ltolerance()} apply to legacy likelihoods and their diagnostic
rows; the fixed public QLE uses its original optimizer settings.

{title:Known reference sign contract}
{p 4 4 2}
Public CNLSDDF returns epsilon=d=alpha+beta'x-gamma'y. Our production residual is
-d. The public generic production StoNED routines use d without converting this
DDF sign to the v-u convention. This can reverse the skewness interpretation.
The reference path preserves that behavior for reproducibility and displays this
qualification. Do not silently label a sign-corrected call as the literal public API.
The historical 2021 frontier mean-shift sign correction is a different change.

{p 4 4 2}
Gamma coefficients can be nonunique even for identical optimal distances. Public
scores use gamma_i'y_i and can therefore depend on which optimal coefficients
were selected. Public parity is conditional on the same coefficient state.
Zero, negative or nonfinite weighted denominators are always refused.
Invalid variance parameters are refused by default. Missing coefficients require an explicit first-stage refit
with {cmd:coefficients(prefix)}. Quantities, directions and saved coefficients
must satisfy the unchanged full-row provenance binding. Sorting is allowed.
No public {cmd:get_stoned()} envelope/target wrapper is supplied.

{title:Raw output policy}
{p 4 4 2}
Default {cmd:invalid(refuse)} rejects inadmissible public parameters and QLE
optimizer nonconvergence. Explicit {cmd:invalid(raw)} retains their raw public
outputs, including negative sigma_u and nonfinite values (Stata missing).
This is an output policy within {cmd:implementation(pystoned)}, not a third
estimator, parameter repair, sign correction or fallback. It is unavailable
with {cmd:implementation(legacy)}. Public negative conditional quantities are
retained under either policy when the parameter/optimizer checks pass.

{p 4 4 2}
Always inspect {cmd:r(parameter_valid)}, {cmd:r(score_valid)} and
{cmd:r(estimate_valid)}. Parameter validity requires finite sigma_u>=0 and
sigma_v>0; score validity requires finite ratios and finite nonnegative raw U.
Combined validity requires both and, for QLE, optimizer success. These checks
do not establish statistical accuracy, the production-residual sign convention,
or uniqueness of the coefficient-dependent ratio. Even when all three flags
equal 1, the raw-d sign and gamma selection qualifications above still apply.
{cmd:r(invalid_reason)} describes
parameter/optimizer problems; {cmd:r(score_status)} separately identifies
negative U or nonfinite scores. Raw numeric return code 0 means the output was
retained successfully, not that the estimate is valid. Text macros preserve
nonfinite parameter tokens. Malformed bridge responses, unavailable Python,
invalid provenance and nonpositive denominators remain errors in raw mode.

{title:Legacy}
{p 4 4 2}
{cmd:implementation(legacy)} retains the stat1 MM/PSL/QLE implementation,
stable directional E(u|z), E(v|z), signed frontier distance and optional targets.
Its outputs are {cmd:prefix_u}, {cmd:prefix_v}, {cmd:prefix_distance}; all original
boundary policies and source checks remain. See the bundled {cmd:stonedddf_legacy.sthlp}
for the complete legacy model, options and results. Its old default syntax now
requires explicit {cmd:implementation(legacy)} through this wrapper.

{p 4 4 2}
Here legacy names the earlier directional implementation; it is not a container
for upstream errors. Its MM uses divisor-N centered moments. The separate scalar
{cmd:stoned, implementation(legacy)} uses a different sample-moment convention;
do not require finite-sample equality merely because both options are named legacy.
The scalar production/VRS {cmd:stoned2012} reference is supplied separately and
does not turn this DDF extension into a reproduction of a 2012 multi-output model.

{title:Diagnostics}
{p 4 4 2}
{cmd:diagnose} returns {cmd:r(diagnose)}, seven rows by twelve columns. Rows:
public_mm, public_qle (literal d); aligned_mm, aligned_qle (counterfactual
-(d-mean(d))); legacy_mm, legacy_psl, legacy_qle (unchanged directional methods).
Aligned rows alter the residual definition and are not original public calls.
They retain the public conditional formula, its 1e-6 denominator adjustment and
the gamma-weighted ratio. Sign alignment alone does not make these quantities
equal to stable legacy conditional inefficiency or to geometric targets.
Legacy score_mean is only a common diagnostic mapping 1-U_stable/(gamma'y),
not a new legacy output or a validated efficiency measure. It is missing if
coefficients or positive denominators are unavailable.

{p 4 4 2}
Columns: rc, sigma_u, sigma_v, mu, U_mean, score_mean, negative_U_n,
denominator_min, mean_distance, converged, reference_sign, parameter_valid.
rc is the Stata return code; failed rows remain explicit. reference_sign=1 for
raw public and -1 otherwise. parameter_valid refers only to nonnegative sigma_u
and positive sigma_v; it is not convergence, score validity or statistical accuracy.
Converged is missing for MOM and a numerical status for likelihoods. On a failed
public QLE row under the default refuse policy, a separate raw probe can retain inadmissible parameters and
optimizer status; its score columns remain missing. Missing Python affects
reference QLE rows only if the selected estimator can run without it.

{p 4 4 2}
The selected {cmd:invalid()} policy also applies to the first four diagnostic
rows. With {cmd:invalid(raw)}, rc=0 can accompany parameter_valid=0.
The original seven-by-twelve {cmd:r(diagnose)} layout is preserved.
New {cmd:r(diagnose_validity)} has four rows (public_mm, public_qle,
aligned_mm, aligned_qle) and three columns (parameter_valid, score_valid,
estimate_valid). Unavailable checks remain missing; a successful raw probe
on a refused QLE row records parameter validity and estimate_valid=0,
while score_valid remains missing. Legacy diagnostic rows retain their original
semantics. {cmd:r(diagnose_invalid_policy)} records the applied policy.

{p 4 4 2}
Diagnostics do not refit stage one, select a method, replace failures, or change
the selected output variables. They require no interpretation as hypothesis tests.
All paths preserve first-stage e(); postestimation results are in r().
A selected public failure does not overwrite pre-existing score variables.

{title:Returned results}
{p 4 4 2}
Public scalars: N, sigma_u, sigma_v, mu, U_mean, score_mean, negative_U_n,
denominator_min, mean_distance, converged. Public macros: cmd, method,
implementation, reference_version, source_status, residual_convention, sign_status,
generate. Legacy returns its original r() plus implementation. With diagnostics:
r(diagnose), r(diagnose_validity), r(diagnose_scope), r(diagnose_invalid_policy).
Public validity scalars: parameter_valid, score_valid, estimate_valid.
Public QLE scalars: qle_optimizer_success, qle_optimizer_status,
qle_optimizer_nfev (missing for MOM). Additional macros: invalid_policy,
invalid_reason, score_status, raw_sigma_u_text, raw_sigma_v_text,
raw_lambda_text, qle_optimizer_message. QLE-only macros are empty for MOM.
No standard errors are returned.

{title:Example}
{p 8 12 2}
{cmd:cnlsddf, inputs(x1 x2) outputs(y1 y2) gx(1 .5) gy(1 2) coefficients(c)}
{p 8 12 2}
{cmd:stonedddf, generate(ref) diagnose pythonexec("/path/to/pinned/python")}
{p 8 12 2}
{cmd:display r(parameter_valid) " " r(score_valid) " " r(estimate_valid)}
{p 8 12 2}
{cmd:matrix list r(diagnose_validity)}
{p 4 4 2}
Read or save these r() results before another command overwrites them.
This is a compatibility example; flags equal to 1 do not resolve the DDF sign
or coefficient-selection issue. For directional quantities, select legacy explicitly:
{p 8 12 2}
{cmd:stonedddf, implementation(legacy) method(qle) generate(own) targets(t)}

{p 8 12 2}
{cmd:stonedddf, method(qle) invalid(raw) generate(raw) diagnose pythonexec("/path/to/pinned/python")}
{p 8 12 2}
{cmd:return list}
{p 8 12 2}
{cmd:matrix list r(diagnose_validity)}
