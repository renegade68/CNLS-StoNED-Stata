{smcl}
{* Public packaging documentation revision 25oct2026; estimator unchanged.}{...}
{* *! version 0.1.0-crs2-dir1 04oct2026}{...}
{title:cnlsddf_sparse — Sparse unregularized joint-input/output CNLS-DDF}

{p 4 4 2}
Optional VRS numerical engine for {help cnlsddf}, plus a separate CRS first-stage
research extension. Both use desirable outputs and a common direction. Requires Stata 16+, macOS/Unix, and a Python environment
with the pinned packages in requirements-tested.txt. Does not change Stata's
Python settings or the original cnlsddf command.


{title:Separate CRS first-stage research candidate}

{p 4 4 2}
This is a research extension beyond the journal scalar commands.
Keep all files in ado/, including cnlsddf_osqp.py and cnlsddf_crs.py, together.
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
CRS support is first-stage only. No CRS decomposition, inference, likelihood optimizer or general speed claim is added.

{title:Syntax}
{p 8 8 2}
{cmd:cnlsddf_sparse} {ifin}, {cmd:inputs(}{it:varlist}{cmd:)}
{cmd:outputs(}{it:varlist}{cmd:)} {cmd:gx(}{it:numlist}{cmd:)}
{cmd:gy(}{it:numlist}{cmd:)} {cmd:python(}{it:executable_path}{cmd:)}
[{cmd:rts(}{it:vrs|crs}{cmd:)} {cmd:distance(}{it:name}{cmd:)} {cmd:residual(}{it:name}{cmd:)}
{cmd:coefficients(}{it:prefix}{cmd:)} {cmd:replace}
{cmd:tolerance(}{it:1e-6}{cmd:)} {cmd:normtolerance(}{it:1e-8}{cmd:)}
{cmd:maxiter(}{it:200000}{cmd:)} {cmd:maxcut(}{it:100}{cmd:)}
{cmd:pivot(}{it:0}{cmd:)} {cmd:full}
{cmd:fullstart(}{it:auto|cold|interior}{cmd:)} {cmd:trace(}{it:path.json}{cmd:)}
{cmd:solver(}{it:osqp|clarabel}{cmd:)} {cmd:ipmmaxiter(}{it:100}{cmd:)}
{cmd:cutbatch(}{it:1}{cmd:)} {cmd:innercap(}{it:1e-3}{cmd:)}
{cmd:innerpolicy(}{it:adaptive|monotone}{cmd:)}
{cmd:workingset(}{it:none|once|bounded}{cmd:)}]

{p 4 4 2}
python() must name an existing executable; use an absolute path. Paths with
quotes, dollar signs, backticks, backslashes or newlines are rejected by the POSIX
shell bridge. Spaces are supported. Dependencies are installed separately.

{p 4 4 2}
Inputs/outputs must be distinct, nonnegative numeric variables. Complete cases
within if/in are used. Directions are common, nonnegative and not all zero,
in the same units/order as quantities. pivot(0) selects the largest scaled
direction component; a positive index must select a positive component.

{p 4 4 2}
Default: solver(osqp) and cutting planes, beginning with adjacent pairs in both directions.
Every scan covers all ordered pairs in blocks, adding at most one new worst
alternative per observation by default. The default workingset(none) does not prune;
all working-set policies keep exhaustive scans.
full instead includes all N(N-1) ordered pairs; it can be faster for small data
but requires quadratic constraint storage. Compare modes on representative data.

{p 4 4 2}
With solver(osqp), maxiter() limits each inner OSQP call; maxcut() limits all calls, including strict
refinement and continuation without new cuts. Both budgets can be reached.
Changed constraints retain only the primal starting point; unchanged QPs reuse
OSQP state. Intermediate tolerances may be relaxed. Final OSQP tolerance is at
most .005*tolerance(); independently recomputed KKT thresholds use .05*tolerance()
with relative scaling. This extra numerical margin is not a universal bound on
fitted-distance error. No objective ridge is used.

{title:CRS first-stage contract}

{p 4 4 2}
{cmd:rts(vrs)} is the default and retains the original VRS calculation.
{cmd:rts(crs)} requires explicit {cmd:solver(clarabel)} and {cmd:workingset(none)}.
Full constraints and cutting planes are supported; OSQP and constraint pruning
are rejected for CRS. Use {cmd:ipmmaxiter()} for the Clarabel iteration budget;
{cmd:maxiter()} applies only to OSQP. AlmostSolved is never certified.

{p 4 4 2}
CRS fixes every alpha to zero. The saved signed distance need not have zero mean;
{cmd:e(mean_distance)} reports its actual mean and {cmd:e(distance_centered)} is 0.
The residual remains negative distance. The solver eliminates direction
normalization and keeps explicit distances with N linking equalities. No ridge
or constant frontier shift is introduced. Coefficients can be nonunique.

{p 4 4 2}
CRS returns {cmd:e(rts)}={cmd:crs} and {cmd:e(source_schema)}={cmd:ddf_crs_rowhash_v1}.
{help stonedddf} rejects CRS with r(301). Changing only e(rts) to vrs does not make
it compatible: the source schema is refused. The VRS constant-shift decomposition
cannot be reused as a CRS model. A signed algebraic distance can project outside
the nonnegative quantity domain and is not automatically a feasible DEA target.

{p 4 4 2}
Additional CRS diagnostics include {cmd:e(max_equality)}, {cmd:e(equalities)},
{cmd:e(free_slopes)} and {cmd:e(rts_crs)}. {cmd:e(parameters)} counts distances plus
free slopes. Sparse matrix bytes include retained setup matrices; the dense matrix
size is a dimensional estimate, not observed resident memory.

{title:Results and postestimation}
{p 4 4 2}
e(converged)=1 requires the selected solver's strict Solved status, finite quantities, the final
inner tolerance, full Afriat feasibility in original distance units, nonnegative
scaled slopes, original-unit direction normalization, and independently checked
primal feasibility, multiplier sign, stationarity and complementarity. OSQP
solved inaccurate and Clarabel AlmostSolved are never final success. Finite noncertified iterates may be
returned with e(converged)=0; {help stonedddf} rejects them. Backend errors leave
previous e() and destination variables intact.

{p 4 4 2}
Default variables are _cnlsddf_distance and _cnlsddf_resid (its negative).
coefficients(c) adds c_alpha, c_beta1,... and c_gamma1,... . Slopes need not be
unique even when fitted distances are unique. Stage-one signed distances are
not inefficiency scores.

{p 4 4 2}
e(cmd)="cnlsddf" maintains the existing postestimation contract;
e(estimator)="cnlsddf_sparse" and e(engine)="sparse_osqp" or "sparse_clarabel"
identify the selected engine. e(qp_solver) is osqp or clarabel.
For certified VRS fits only, use {cmd:stonedddf, implementation(legacy) method(mm)}, {cmd:method(psl)} or {cmd:method(qle)} immediately
after a certified fit, subject to the same distributional and boundary checks.

{p 4 4 2}
Saved scalars include N, m, s, rss, mean_distance, converged, pivot, tolerance,
normtolerance, max_viol, max_mono, max_mono_scaled, max_normalization, primal,
dual_feasibility, stationarity, complementarity and corresponding eps_* thresholds,
solver_status, inner_tolerance, total_inner_iters, cp_iters, n_active, setups,
same_qp_resumes, expanded_primal_starts, parameters and objective_ridge (zero).
e(gx), e(gy), e(cp_history) and e(sample) are also saved.

{p 4 4 2}
e(seconds) times the Python fit, excluding Python startup and file transfer.
A Stata timer around the command includes that overhead. matrix_payload_bytes
counts the CSC arrays of Z, P and A; dense_normal_matrix_bytes is the hypothetical
8(N(m+s))^2 dense normal matrix avoided by this implementation. Neither is total
or peak process memory; solver factorization fill-in is not counted.
scan_workspace_bytes describes one block-by-N double array, not all temporaries.

{title:Example}
{phang2}{cmd:cnlsddf_sparse, inputs(x1 x2) outputs(y1 y2) gx(1 0) gy(0 1) python("/absolute/path/.venv/bin/python") coefficients(c) full}{p_end}
{phang2}{cmd:stonedddf, implementation(legacy) method(psl) generate(sd)}{p_end}
{phang2}{cmd:cnlsddf_sparse, inputs(x1 x2) outputs(y1 y2) gx(1 0) gy(0 1) python("/absolute/path/.venv/bin/python") solver(clarabel) ipmmaxiter(100) maxcut(100) replace}{p_end}

{title:Direct interior-point engine (dev3)}
{p 4 4 2}
{cmd:solver(clarabel)} solves each restricted QP in exact distance coordinates
using Clarabel 0.11.1, then maps the solution and multipliers back to the original
QP for the same independent certificate. The objective has zero ridge. Each scan
still covers all ordered Afriat pairs; the default adds one worst new
alternative per observation. Pruning is optional via workingset(); there is no
warm starting or factor reuse: every inner call constructs a fresh interior-point model.

{p 4 4 2}
{cmd:ipmmaxiter()} limits iterations per direct Clarabel call; {cmd:maxiter()} is
an OSQP control and has no effect on this engine. {cmd:maxcut()} bounds all inner
calls, including strict refinement without new constraints. An inner iteration
limit returns an uncertified result; there is no automatic fallback to OSQP.
Finite AlmostSolved iterates may inform intermediate cuts, but can never certify
the final estimator. Other unsuccessful native states stop with diagnostics,
subject to the one-shot restoration rule under optional workingset() pruning.

{p 4 4 2}
Intermediate tolerances may be relaxed after cuts are added. Final feasibility
and gap tolerances are at most .005*tolerance(), and strict Solved status plus
the unchanged original-coordinate KKT/shape checks are mandatory. A Solved iterate
that fails the external checks can be solved again at a tighter tolerance within
maxcut(). Merely finding no new violated pair does not certify success.

{p 4 4 2}
{cmd:solver(clarabel) full} directly solves the same full QP and applies the same
certificate. It does not run an OSQP refinement stage. Explicit fullstart(cold)
or fullstart(interior) is invalid for this engine; fullstart(auto) is accepted as
the default and e(full_start) reports none. e(initializer) is none and initializer
scalars are zero: Clarabel is the selected solver, not an initializer here.

{p 4 4 2}
With this engine, e(total_inner_iters) and e(ipm_iterations) count interior-point
iterations; e(osqp_iterations)=0. e(ipm_maxiter), e(constraint_expansions),
e(strict_refinements), e(ipm_matrix_bytes) are additionally saved.
same_qp_resumes, expanded_primal_starts and status_resets are zero.
ipm_matrix_bytes is the peak T/Q/G array payload, not peak process memory.
The seven-column cp_history format is unchanged, but solver_status is the native
status of the selected solver; consult e(qp_solver) when interpreting it.
Trace JSON identifies the solver and contains per-QP native residuals, objectives,
iteration records and logs, alongside the independent certificate diagnostics.

{title:Restricted-QP policy controls (dev4/dev5)}
{p 4 4 2}
For {cmd:solver(clarabel)} cutting planes, {cmd:cutbatch()} sets the maximum number
of newly violated alternatives retained per observation per scan. Candidates are
ordered by descending violation, with the lowest index breaking ties. Active
pairs and self-pairs are excluded. The maximum violation and total violating-pair
count are computed before candidate filtering and remain exhaustive for every
batch size. More cuts can reduce QP calls but increase storage and factorization.

{p 4 4 2}
With the default {cmd:innerpolicy(adaptive)}, {cmd:innercap()} sets the intermediate feasibility/gap tolerance ceiling after
adding cuts. The requested tolerance is max(.005*tolerance(),
min(innercap(), .01*global_violation)). It does not relax the final certificate.
A smaller ceiling can improve cut selection at the cost of more inner iterations.
The first QP is strict; final strict Solved and original-coordinate KKT/global
shape checks are unchanged. These tolerance controls add no warm starts or
objective penalty. Pruning is controlled separately by workingset().

{p 4 4 2}
{cmd:innerpolicy(monotone)} permits the initial strict QP to release its tolerance
when adding the first cuts. From QP 2 onward, the next requested tolerance cannot
exceed the current tolerance, including after further cuts. Thus all subsequent
tightenings are retained. The rule applies from QP 2, not only after reaching the
final strict tolerance. It resets on each estimate, and does not reuse an IPM
iterate. Refinements below the final tolerance are also retained. This option is
available only for Clarabel CP; it can increase inner work and is not a universal
speed improvement. Default {cmd:innerpolicy(adaptive)} reproduces dev4.
e(inner_policy) and the JSON top-level inner_policy record the policy; for Clarabel,
e(inner_monotone) is its 0/1 indicator and e(tolerance_holds) counts prevented
tolerance increases after adding cuts. The seven-column history is unchanged.

{p 4 4 2}
Defaults {cmd:cutbatch(1) innercap(1e-3)} reproduce dev3. A combined option to
compare on representative data is {cmd:cutbatch(4) innercap(1e-5)}.
Nondefault settings require Clarabel CP; they are rejected with OSQP or full.
cutbatch() must be a positive integer; innercap() must be in (0,1e-3]. The strict
final tolerance is a floor even if innercap() is smaller. Python adaptive=False
continues to request the strict tolerance regardless of the ceiling.
e(cut_batch) and e(inner_cap) record settings for the direct Clarabel engine.
The candidate count in cp_history can now exceed N; its seven-column layout is
unchanged. maxcut() still bounds QP calls, not the number of individual pairs.

{phang2}{cmd:cnlsddf_sparse, inputs(x1 x2) outputs(y1 y2) gx(1 1) gy(1 1) python("/absolute/path/.venv/bin/python") solver(clarabel) cutbatch(4) innercap(1e-5)}{p_end}

{title:Optional working-set pruning and recovery (dev6)}
{p 4 4 2}
{cmd:workingset(none)} is the default and preserves the dev5 solver trajectories.
{cmd:workingset(once)} or {cmd:workingset(bounded)} requires
{cmd:solver(clarabel)} without {cmd:full}. They can be combined with cutbatch(),
innercap(), and either innerpolicy(). Pruning is a numerical option for the same
unregularized model, not a statistical restriction on which Afriat pairs count.

{p 4 4 2}
once permits removal after QP 10; bounded permits it after QPs 5, 10, 15, and 20.
Removal occurs only when new cuts are being added, native status is Solved,
the requested inner tolerance is looser than the final threshold, and another
QP fits within maxcut(). Initial adjacent pairs are protected. At most half of
the active pairs are removed per event. A pair can be removed at most once;
violated omitted pairs can re-enter through the exhaustive scan and then stay.

{p 4 4 2}
Selection uses original distance-row slack and multipliers. Slack must exceed
10*max(tolerance,inner_tolerance)*(1+max(abs(own_distance),abs(alternative_distance))),
and abs(multiplier) must not exceed
max(tolerance,inner_tolerance)*(1+max(abs(distance))). Eligible pairs are ranked
by descending slack, then row index. These fixed heuristics can help some data
and worsen other cases; compare time, memory, and certification on your data.

{p 4 4 2}
Pruning enables at most one restoration attempt after a finite native
NumericalError/InsufficientProgress, or a terminal AlmostSolved with no new
cut candidates at the final inner tolerance. The attempt requires currently
omitted previously removed pairs and room within maxcut(). Existing row order
is retained and missing removed pairs are appended in lexicographic order.
Subsequent pruning is disabled and all later QPs use at least the final accuracy.
New violated pairs may still be added. Restoration does not add all N(N-1)
constraints, does not reset maxcut()/ipmmaxiter(), and does not guarantee success.
MaxIterations, MaxTime, infeasibility, and exceptions do not trigger recovery.
Solved results that fail external certification retain the existing refinement
path. AlmostSolved never counts as final success, even after restoration.

{p 4 4 2}
e(working_set) records none/once/bounded. New scalars are working_set_enabled,
working_set_prunes (removal events), working_set_removed (total removed pairs),
working_set_readmitted (violation-driven returning pairs), working_set_pinned,
recovery_enabled, recovery_attempts (0 or 1), and recovery_restored.
Restoration counts are separate from violation-driven readmission/pinning.
All new counters/indicators are zero under workingset(none), including OSQP/full.
A finite uncertified result still posts e(converged)=0 and cannot feed stonedddf;
backend/option/trace errors preserve previous e() and destination variables.

{p 4 4 2}
The seven-column e(cp_history) remains unchanged, and includes failed/recovery QPs.
Trace JSON adds working_set, working_set_history (ordered removals/additions),
and recovery_history (trigger, source iteration, restored pairs, next row count,
and requested tolerance). Detailed pair ledgers are collected only with trace();
without it the Python result retains recovery summaries but no pair lists.
The default none policy keeps pair ledgers empty even with trace().

{phang2}{cmd:cnlsddf_sparse, inputs(x1 x2) outputs(y1 y2) gx(1 1) gy(1 1) python("/absolute/path/.venv/bin/python") solver(clarabel) cutbatch(4) workingset(bounded) trace("working_set.json")}{p_end}

{title:Full-constraint initialization and traces (dev2)}
{p 4 4 2}
With {cmd:solver(osqp) full}, {cmd:fullstart(auto)} initializes from Clarabel 0.11.1 solving
all constraints in an exactly equivalent distance-coordinate QP. This is a full
interior-point solve, not a cutting-plane fit. It adds no objective ridge.
The final OSQP solve still includes every ordered pair and must pass the same
strict status, tolerance and original-coordinate certificate. The initializer
alone never declares success.

{p 4 4 2}
{cmd:fullstart(cold)} skips initialization for comparison. {cmd:fullstart(interior)}
requires {cmd:full}. Auto without full keeps the existing CP start. Clarabel is
loaded only when needed; dependency/setup errors are reported. A finite initializer
that lacks strict Solved status or fails the original certificate is discarded,
then the full OSQP solve starts cold. Its feasibility and gap tolerances are .00005*tolerance(), 100 times stricter
than the final OSQP setting, to provide a sufficiently accurate start. The final
original-coordinate certification thresholds remain unchanged. Initializer iterations are capped at
min(maxiter(),100), separately from the OSQP maxiter()/maxcut() budgets.

{p 4 4 2}
{cmd:trace("path.json")} writes solver diagnostics and native iteration logs.
An existing file requires {cmd:replace}. The JSON contains Clarabel iteration
records, per-call OSQP residuals, rho update counts/estimates, certificate metrics,
and native OSQP lines with the current rho (normally every 200 iterations and
when rho changes). rho_estimate is a proposed value, not necessarily current rho.
Trace-file failures occur before estimation results or destination variables change.

{p 4 4 2}
New e() scalars: initializer_used, initializer_status (Clarabel raw status),
initializer_iterations, initializer_seconds, initializer_matrix_bytes,
initializer_tolerance, status_resets.
Initializer status 1 means Solved; status 0 with zero time/iterations means not
attempted. e(initializer) identifies clarabel/none; e(full_start) reports the
requested effective policy, even if an unsuccessful initializer was discarded.
e(trace_file) contains the requested path. Under solver(osqp), total_inner_iters counts OSQP only;
initializer_iterations counts interior-point steps separately. seconds includes
both phases. initializer_matrix_bytes counts T/Q/G arrays during initialization,
not peak memory. Existing seven-column e(cp_history) remains unchanged.

{p 4 4 2}
Each continuation updates the identical linear objective vector to reset OSQP's
information fields. In the pinned version, settings-only continuation can retain
an old solved status. The q-only update preserves the mathematical problem and
warm iterate/factors; independent certification remains mandatory.

{title:CRS example}

{p 4 4 2}
With nonnegative x1 x2 y1 y2 already loaded and a direction chosen in their units:

{phang2}{cmd:cnlsddf_sparse, inputs(x1 x2) outputs(y1 y2) gx(1 1) gy(1 1) rts(crs) solver(clarabel) python(/absolute/path/.venv/bin/python) coefficients(c)}{p_end}
{phang2}{cmd:assert e(converged)==1}{p_end}
{phang2}{cmd:assert c_alpha==0 if e(sample)}{p_end}

{title:Direction precision}

{pstd}
Explicit numeric lists in {cmd:gx()} and {cmd:gy()} preserve their parsed double
precision and are sent to Python with 17 significant digits. The saved
{cmd:e(gx)} and {cmd:e(gy)} use these same values. If a direction uses compact
range notation, such as {cmd:1/3} or {cmd:0(1)2}, that entire direction retains
Stata's numlist expansion and its precision. Use an explicit list when the
last digits matter; expand ranges before calling the command if necessary.

{pstd}
{cmd:e(gx_syntax)} and {cmd:e(gy_syntax)} report {cmd:literal} or {cmd:numlist};
{cmd:e(bridge_version)} is {cmd:0.1.0-dir1}. The original dir1-only candidate retained the crs1 engine; this candidate uses
the crs2 numerical profile described below.
This interface correction does not change the optimization tolerance or make
an uncertified CRS solution converged. CRS postestimation remains blocked.

{title:CRS numerical profile in this research candidate}

{pstd}
The crs2 backend uses internal KKT static regularization 1e-6 and iterative
refinement absolute/relative tolerances 1e-14 with at most 30 refinement steps.
This does not add a ridge to the CNLS objective. The existing strict native and
original-coordinate certificate, exhaustive Afriat scan, and constraint-addition
loop remain in force. VRS numerical settings are unchanged. CRS reports
{cmd:e(kkt_profile)} as {cmd:crs2_reg6_refine}; JSON diagnostics record the settings.
The Stata direction parser remains dir1. This isolated research candidate is not
a production release and does not enable CRS StoNED postestimation.
