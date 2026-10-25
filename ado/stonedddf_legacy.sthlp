{smcl}
{* Public packaging documentation revision 25oct2026; estimator unchanged.}{...}
{* *! help revision ddf-interpretation1 25oct2026; executable stat1 unchanged}{...}
{title:stonedddf -- Directional StoNED MM, PSL, and QLE postestimation}

{p 4 4 2}
Experimental normal/half-normal MM, PSL, or QLE decomposition after certified VRS {help cnlsddf}.
Supports joint multiple inputs and desirable outputs with a common direction.
Preserves the stage-one {cmd:e()} results and returns stage-two results in {cmd:r()}.


{title:Separate safety research candidate}

{p 4 4 2}
This is a research extension beyond the journal scalar commands.
Keep all files in ado/ together.
Stage one records {cmd:e(source_schema)}, {cmd:e(source_signature)}, and
{cmd:e(source_variables)}. Stage two requires this binding; older fits must be refitted.
It checks the complete in-sample rows (inputs, outputs, distance, residual and any
saved coefficients) together with the direction. Sorting, display formats and
out-of-sample edits are allowed. Changed or missing source data are refused with r(459).
The checksum detects ordinary changes; it is not a cryptographic certificate.

{p 4 4 2}
Stage-two output replacement uses exact variable names and restores the dataset
on late write errors. Successful calls return {cmd:r(source_status)} =
{cmd:ddf_rows_verified} and {cmd:r(source_schema)} = {cmd:ddf_rowhash_v1}.
This does not add CRS, inference, a new likelihood optimizer or a general speed claim.

{title:Syntax}

{p 8 12 2}
{cmd:stonedddf}, {cmd:implementation(legacy)} [ {opt gen:erate(prefix)} {opt tar:gets(prefix)}
{opt noisea:djusted(prefix)} {opt replace}
{opt method(mm|psl|qle)} {opt maxiter(#)} {opt ltol:erance(#)}]

{p 4 4 2}
Run immediately after {cmd:cnlsddf}, using its entire {cmd:e(sample)} (at least
three observations). No varlist, weights, {cmd:if}, or {cmd:in} are accepted.
The default method is {cmd:mm}. The default {cmd:generate()} prefix is {cmd:_stonedddf}.
All outputs are double precision and missing outside the estimation sample.

{title:Outputs and interpretation}

{p 4 4 2}
{cmd:generate(p)} creates {cmd:p_u} = E(u|z), {cmd:p_v} = E(v|z), and
{cmd:p_distance} = signed distance D to the shifted frontier, where z = -D = v-u.
Inefficiency {cmd:p_u} is nonnegative. Distance D and noise {cmd:p_v} are signed.
These are in directional distance units; none is a Farrell ratio or an
exponentiated efficiency score. Changing the magnitude of the direction changes
these units.

{p 4 4 2}
{cmd:targets(t)} creates {cmd:t_x1}, ..., {cmd:t_y1}, ... in the input/output
order of {cmd:cnlsddf}. The target is (x-D*gx, y+D*gy). It lies on the shifted
piecewise-linear frontier up to the stage-one numerical tolerances. A negative
coordinate means that this projected point falls outside the nonnegative
physical domain; {cmd:t_negative}=1 flags that observation. Values are not clipped.

{p 4 4 2}
{cmd:noiseadjusted(n)} creates {cmd:n_x1}, ..., {cmd:n_y1}, ... for
(x+E(v|z)*gx, y-E(v|z)*gy), and {cmd:n_negative} flags any negative coordinate.
Moving this noise-adjusted point by E(u|z) along (-gx,+gy) gives the target.
Moving the observed point by E(u|z) alone generally does not reach the frontier.

{p 4 4 2}
{cmd:replace} permits replacement of existing output variables only after
successful calculation. Output names must be distinct and cannot overwrite
stage-one input, output, distance, residual, or recorded coefficient variables.
Stage-one coefficients need not have been saved to calculate targets.

{title:Model, moments, and numerical behavior}

{p 4 4 2}
Assume independent v~N(0,sigma_v^2), u~|N(0,sigma_u^2)| with additive error
along the supplied common direction. This is a scalar directional error model,
not unrestricted multivariate noise. Let d be the stage-one signed distance,
dbar its sample mean, and r=-(d-dbar). Both central moments use denominator N:
m2=mean(r^2), m3=mean(r^3). For negative m3,
sigma_u=(-m3/[(4/pi-1)*sqrt(2/pi)])^(1/3) and
sigma_v^2=m2-sigma_u^2*(1-2/pi).

{p 4 4 2}
The label legacy selects this earlier directional estimator, not erroneous
upstream code. Scalar {cmd:stoned, implementation(legacy)} MM instead uses a
sample-variance convention. That finite-N difference is not a DDF sign error.
For a single output with gx=0 and gy=1, the divisor-N MM definition and
observed-point frontier shift agree with the separate scalar {cmd:stoned2012}
MOM reference under the same centered production/VRS fit. This limited reduction
does not validate general multiple-output inference or assert that the DDF
likelihood alternatives reproduce the 2012 concentrated likelihood.

{p 4 4 2}
Let mu=sigma_u*sqrt(2/pi). The frontier shifts by mu-dbar, D=d-dbar+mu,
and z=-D. E(u|z) is the conditional truncated-normal mean; E(v|z)=z+E(u|z).
Extreme negative standardized conditional means use a continued fraction to
avoid normal-tail underflow and cancellation. Missing or nonfinite results
cause an error instead of being converted to apparent efficiency.

{p 4 4 2}
For method(mm), positive m3 has the wrong production skew. The explicit boundary fallback is
sigma_u=0, sigma_v=sqrt(m2), status {cmd:wrong_skew}; this does not match the
third moment or establish universal efficiency. Exactly zero m3 uses
{cmd:zero_skew}. Exactly zero variance is rejected with r(498); it does not identify a supported decomposition.
If negative m3 implies nonpositive noise variance, the command stops with
r(498), without replacing outputs. The positive-inefficiency zero-noise limit is not implemented.
The likelihood methods described below offer alternatives when MM fails.

{p 4 4 2}
Do not modify in-sample estimation data between stages. In addition to sample
size, missingness, residual sign, RSS and distance mean, the safety1 protection
checks complete input/output/distance/residual rows, any saved coefficients and
the common direction. It allows sorting, unrelated columns and out-of-sample edits.
This ordinary-change checksum is not cryptographic authentication. Refit after
changing stage-one data. Numerical convergence does not establish statistical
identification or validity of distributional assumptions. No standard errors are provided.

{title:Likelihood alternatives}

{p 4 4 2}
{cmd:method(psl)} estimates both scales by maximizing the centered
normal/half-normal pseudo-likelihood. {cmd:method(qle)} maximizes the same
criterion subject to sigma_v^2 + (1-2/pi)*sigma_u^2 = m2, reducing it to
one parameter. QLE is a moment-constrained curve, not the exact scale profile
of the two-parameter PSL objective. Both use z=r-sigma_u*sqrt(2/pi), with r
centered as above. They are distinct estimators and need not agree with MM
or with each other. They can be used even when MM implies negative noise variance.

{p 4 4 2}
PSL profiles the reciprocal total scale using a strictly convex one-dimensional
problem at each variance ratio. The outer angle atan(sigma_u/sigma_v) is searched
on a 65-point grid, with interval refinement of detected minima and both endpoint
neighborhoods. QLE uses the same outer search with its scale fixed by m2.
{cmd:ltolerance()} defaults to 1e-8: it controls angle interval width and relative
inner scale interval width. Boundary likelihood comparisons use a separate fixed
per-observation tolerance of 1e-8, regardless of ltolerance().
It must be positive and no greater than .001. {cmd:maxiter()} defaults to 200
per scalar refinement (not a total evaluation budget). These controls do not
change MM calculations. Search completion is numerical, not a proof of a global
optimum for every possible sample. Inspect {cmd:r(likelihood_profile)} when needed.

{p 4 4 2}
The exact sigma_u=0 normal model is an explicit candidate. If its mean criterion
is within the fixed 1e-8 mean-criterion tolerance of the best candidate, {cmd:normal_boundary} is
returned with {cmd:r(boundary)}=4. This is not evidence of universal efficiency.
No skew-sign shortcut forces this outcome. Likelihood scales with zero sample
variance have no finite optimum and are rejected. MM also rejects zero variance.

{p 4 4 2}
The positive-noise search stops at sigma_v/sqrt(sigma_u^2+sigma_v^2)=1e-5.
If that endpoint is competitive with an improving interior fit at the
fixed 1e-8 mean-criterion tolerance, or its angle is within 4*ltolerance() of
the upper search endpoint, r(498) is issued. This may indicate a
zero-noise support-boundary optimum or an insufficiently separated near-boundary
fit. This version does not certify that nonregular model. The noise floor is
not silently returned as an estimate. Iteration exhaustion issues r(430).
Neither failure falls back silently to MM or changes requested output variables.

{title:Returned results}

{p 4 4 2}
Scalars: {cmd:r(N)}, {cmd:r(sigma_u)}, {cmd:r(sigma_v)}, {cmd:r(mu)},
{cmd:r(mean_distance)}, {cmd:r(frontier_shift)}, {cmd:r(m2)}, {cmd:r(m3)},
{cmd:r(boundary)}, {cmd:r(target_negative)}, {cmd:r(noise_negative)}.
The last two count flagged observations, or are missing when the respective
outputs were not requested. Boundary codes: 0 interior, 1 wrong skew,
2 zero skew, 4 likelihood normal boundary. Code 3 is reserved for historical results;
this candidate rejects zero variance instead of returning it.

{p 4 4 2}
Macros: {cmd:r(cmd)} = stonedddf, {cmd:r(method)} = mm, psl, or qle,
{cmd:r(distribution)} = normal_halfnormal, {cmd:r(status)},
{cmd:r(generate)}, {cmd:r(targets)}, {cmd:r(noiseadjusted)}.
As with other r-class commands, subsequent commands may replace {cmd:r()}.
The unchanged {cmd:e(cmd)} remains {cmd:cnlsddf}.

{p 4 4 2}
Likelihood methods additionally return {cmd:r(ll)} (including density constants
and original distance units), {cmd:r(converged)}=1 on success,
{cmd:r(search_iterations)} (outer interval contractions), {cmd:r(evaluations)}
(total objective/scale-score evaluations), {cmd:r(noise_ratio)},
{cmd:r(noise_floor)}=1e-5 and {cmd:r(ltolerance)}.
The additional {cmd:r(selection)} is {cmd:interior} or {cmd:normal_boundary}.
{cmd:r(normal_ll)} reports the exact normal candidate;
{cmd:r(candidate_ll)} reports the best grid/refined candidate before the tie rule
(including the normal candidate). {cmd:r(ll_gain)} equals their difference divided
by N. {cmd:r(boundary_tolerance)} is the fixed numerical tie threshold 1e-8,
not a significance test or an identification threshold. {cmd:r(ll)} always belongs
to the selected parameters. These likelihood diagnostics are not retained after MM.
{cmd:r(likelihood_profile)} is a 65-by-4 matrix with columns
{cmd:tau sigma_u sigma_v ll}, where tau=sigma_u/sqrt(sigma_u^2+sigma_v^2).
Reported convergence means the prescribed grid/refinement and boundary checks
completed. After any error, no new stage-two results are certified.

{title:Example}

{phang2}{cmd:cnlsddf, inputs(x1 x2) outputs(y1 y2) gx(1 0) gy(0 1)}{p_end}
{phang2}{cmd:stonedddf, implementation(legacy) generate(sd) targets(target) noiseadjusted(clean)}{p_end}
{phang2}{cmd:return list}{p_end}
{phang2}{cmd:stonedddf, implementation(legacy) method(qle) generate(qle) targets(qtarget)}{p_end}
{phang2}{cmd:stonedddf, implementation(legacy) method(psl) generate(psl) targets(ptarget)}{p_end}

{p 4 4 2}
Choose the direction in the original quantity units before fitting. An example
does not ensure that a particular dataset has admissible MM variance estimates.

{title:References and scope}

{p 4 4 2}
Dai, Fang, Lee, and Kuosmanen (2024), pyStoNED: A Python Package for Convex
Regression and Frontier Estimation, Journal of Statistical Software 111(6).
{browse "https://doi.org/10.18637/jss.v111.i06":doi:10.18637/jss.v111.i06}.
This is a separate implementation, not a claim of bitwise compatibility
with pyStoNED or scalar {help stoned_legacy}. See {help cnlsddf_applied} for
direction and target interpretation. Panel extensions, bootstrap inference,
alternative noise models and Simar-Wilson estimation are not supplied.
