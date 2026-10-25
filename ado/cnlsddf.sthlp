{smcl}
{* Public packaging documentation revision 25oct2026; estimator unchanged.}{...}
{* *! version 0.1.0-safety1 01oct2026}{...}
{title:cnlsddf — Joint-input/output convex nonparametric directional distance}

{p 4 4 2}
Stage-one CNLS with variable returns to scale, desirable outputs and a common
nonnegative direction. Requires Stata 16. Experimental development implementation.


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
{cmd:cnlsddf} {ifin}, {cmd:inputs(}{varlist}{cmd:)}
{cmd:outputs(}{varlist}{cmd:)} {cmd:gx(}{it:numlist}{cmd:)} {cmd:gy(}{it:numlist}{cmd:)}
[{cmd:distance(}{it:name}{cmd:)} {cmd:residual(}{it:name}{cmd:)}
{cmd:coefficients(}{it:prefix}{cmd:)} {cmd:replace}
{cmd:tolerance(}{it:#}{cmd:)} {cmd:normtolerance(}{it:#}{cmd:)}
{cmd:maxiter(}{it:#}{cmd:)} {cmd:maxcut(}{it:#}{cmd:)}
{cmd:rho(}{it:#}{cmd:)} {cmd:pivot(}{it:#}{cmd:)} {cmd:full} {cmd:cpstate(}{it:cold|mapped}{cmd:)} {cmd:innerpolicy(}{it:adaptive|strict}{cmd:)} {cmd:verbose}]

{title:Model and interpretation}

{p 4 4 2}
Minimize one half of the sum of squared signed distances
{it:d_i = alpha_i + beta_i'x_i - gamma_i'y_i}, subject to
{it:d_i <= alpha_h + beta_h'x_i - gamma_h'y_i} for every pair,
{it:beta_i'gx + gamma_i'gy = 1}, and nonnegative beta and gamma.
The fitted joint surface is the minimum of these affine functions.

{p 4 4 2}
The saved distance is centered stage-one distance; it can be negative.
The saved production residual is {it:r=-d}. With one output, gx=0 and gy=1,
this residual matches production VRS {help convnls}.
Neither variable is a nonnegative inefficiency estimate. This command does not
implement StoNED distributional decomposition, efficiency scores or inference.
Do not automatically interpret exp(-d) or exp(-r) as efficiency.

{title:Data and directions}

{p 4 4 2}
All input/output variables must be distinct and numeric. Complete cases from
{cmd:if}/{cmd:in} and every input/output variable define {cmd:e(sample)}.
At least two complete cases are required. Negative quantities in this sample
are rejected; missing quantities are excluded. Zero and constant columns are allowed.

{p 4 4 2}
Directions are required, finite, nonnegative and not jointly zero. Their lengths
and order must match the corresponding variable lists. Input-only, output-only
and sparse directions are supported. Direction is (decrease inputs, increase outputs).
When changing quantity units, rescale the corresponding direction components too.
Multiplying the entire direction by c>0 divides distances by c and RSS by c^2.
CRS, undesirable outputs, observation-specific directions, weights, contextual
variables, cost orientation and logarithmic error models are not supported.

{title:Options}

{p 4 8 2}
{cmd:distance()} saves signed distance; default {cmd:_cnlsddf_distance}.
{cmd:residual()} saves its negative; default {cmd:_cnlsddf_resid}.
Generated variables are double precision and missing outside the estimation sample.

{p 4 8 2}
{cmd:coefficients(prefix)} optionally saves {it:prefix}{cmd:_alpha},
{it:prefix}{cmd:_beta1} ... and {it:prefix}{cmd:_gamma1} ... in original units.
Coefficients can be nonunique even when fitted distances are unique.
Their minimum defines an envelope for prediction, but no predict command is supplied
in this version. Off-sample predictions can depend on the selected optimal coefficients.

{p 4 8 2}
{cmd:replace} permits overwriting named destination variables. Without it, a
collision is an error. Input/output variables can never be overwritten. Destination
names must be distinct and at most 32 characters, including coefficient suffixes.
Validation and estimation run before any existing destination is replaced.
An uncertified finite final iterate is still saved; inspect {cmd:e(converged)}.

{p 4 8 2}
{cmd:tolerance()} defaults to 1e-6. It bounds the exhaustive Afriat violation in
original distance units and the nonnegativity violation in internally scaled slope
units. Internal data scaling divides each quantity column and its direction by the
column maximum absolute quantity (1 for an all-zero column), without centering.
The initial restricted-QP tolerance is 0.05 times the requested tolerance; it may
be tightened if the QP passes but the global shape checks fail.

{p 4 8 2}
{cmd:normtolerance()} defaults to 1e-8 and separately bounds the dimensionless
normalization error after recovering coefficients in original units. No relaxed
30-times-tolerance success state is used. Both tolerances must be positive and finite.

{p 4 8 2}
{cmd:maxiter()} defaults to 20000 iterations per restricted QP call.
{cmd:maxcut()} defaults to 100 solved-and-scanned outer iterations. The total inner
budget can therefore be their product. A budget exit never implies convergence.

{p 4 8 2}
{cmd:rho()} defaults to 1 and must be positive and finite. It initializes the
ADMM penalty; normalized residual balancing subsequently adapts it.

{p 4 8 2}
{cmd:pivot()} defaults to 0 (choose the largest internally scaled direction
component). A positive index into the concatenated gx,gy vector must select a
positive direction component. This eliminates one slope using the normalization
identity and preserves its nonnegativity via an affine upper bound. Alternative
pivots are primarily useful for numerical validation and can alter speed.

{p 4 8 2}
{cmd:full} starts with all N(N-1) Afriat pairs, useful for small reference problems.
The default starts with adjacent bidirectional pairs in observation order and adds
at most one new worst violation per observation after each exhaustive scan.
No constraint pruning is performed. {cmd:verbose} displays iteration diagnostics.

{p 4 8 2}
{cmd:cpstate(cold)} is the default: reset slack and dual variables when pairs
are appended. With {cmd:innerpolicy(strict)}, this matches dev1.
{cmd:cpstate(mapped)} preserves existing rows' slack and multipliers,
initializes new rows with zero dual and projected slack, and resets the penalty
adaptation schedule for the enlarged QP. It accepts only unchanged scaling, bounds
and old constraint order. The objective and returned primal iterate must be unchanged.
{cmd:coldstart} is an alias for cold and cannot be combined with mapped.
Same-QP continuation remains enabled in both modes; {cmd:full} has no pair additions.
Transfer is an initialization, never a convergence certificate. All checks are rerun.

{p 4 8 2}
{cmd:innerpolicy(adaptive)} is the default. The initial solve uses 0.05*tolerance.
After adding pairs, the next intermediate QP tolerance is
max(0.05*tolerance, min(1e-3, 0.01*previous_global_violation)). When no new pair is
found, it returns immediately to the strict tolerance and resumes the same QP.
A successful result always requires an inner tolerance at most 0.05*tolerance,
plus every original convergence check. {cmd:innerpolicy(strict)} keeps all
restricted solves strict, reproducing dev1 with {cmd:cpstate(cold)}. The full mode
has no additions and therefore receives no intermediate tolerance relaxation.

{title:Convergence and saved results}

{p 4 4 2}
{cmd:e(converged)=1} requires original-objective QP stationarity, primal/dual
residuals and inequality complementarity, plus a complete original-unit Afriat scan,
scaled monotonicity and normalization checks. {cmd:e(converged)=0} is uncertified.
{cmd:e(stop_reason)} is certified, maxcut, inner_maxiter, normalization or
numerical_failure. These are numerical tolerance checks, not inferential guarantees.

{p 4 4 2}
Saved scalars include {cmd:N m s rss mean_distance resid_sign pivot parameters},
{cmd:cp_iters n_active total_inner_iters inner_solves warm_starts expanded_starts final_rho},
{cmd:factorizations chol_solves lu_fallbacks rho_updates},
{cmd:max_viol max_mono max_mono_scaled max_normalization tolerance normtolerance},
and {cmd:inner_status inner_iter inner_tolerance}. The four inner diagnostics
{cmd:inner_primal inner_dual inner_stationarity inner_complementarity} have
corresponding {cmd:inner_eps_*} thresholds. {cmd:max_mono} is in original coefficient
units; {cmd:max_mono_scaled} is used for certification.

{p 4 4 2}
Matrices {cmd:e(gx)} and {cmd:e(gy)} store original directions with variable names.
{cmd:e(cp_history)} records iteration, active pairs, inner status, global maximum
violation, candidate count, inner iteration count and inner tolerance. Candidates
on a terminal row are not added.
Macros include {cmd:cmd inputs outputs rts engine stop_reason distance_var resid_var
coefficient_prefix expansion inner_policy}.
{cmd:expanded_starts} counts accepted transfers after adding pairs;
{cmd:warm_starts} counts all subsequent state reuse. {cmd:expansion} is mapped or cold.

{p 4 4 2}
{cmd:constraint_bytes} measures compact operator array payload only;
{cmd:dense_constraint_bytes} is its hypothetical dense G equivalent.
{cmd:normal_matrix_bytes} is one dense P by P matrix, where P=N(m+s).
These are not process memory measurements. Several dense normal/factor matrices,
the pair-membership matrix and temporary arrays are also allocated. Large-sample
scalability is not established by the small reference tests.

{title:Example}

{p 8 8 2}{cmd:cnlsddf, inputs(x1 x2) outputs(y1 y2) gx(1 0) gy(0 1) coefficients(ddf)}
{p 8 8 2}{cmd:display e(converged)}
{p 8 8 2}{cmd:matrix list e(gx)}

{title:References}

{p 4 4 2}
Dai, Fang, Lee and Kuosmanen (2024), pyStoNED: A Python Package for Convex
Regression and Frontier Estimation, Journal of Statistical Software 111(6),
section 4.4. DOI: 10.18637/jss.v111.i06.

{p 4 4 2}
Kuosmanen and Johnson (2017), Modeling Joint Production of Multiple Outputs in
StoNED: Directional Distance Function Approach, European Journal of Operational
Research 262(2), 792-801. DOI: 10.1016/j.ejor.2017.04.014.
