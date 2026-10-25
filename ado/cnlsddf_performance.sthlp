{smcl}
{* Public packaging documentation revision 25oct2026; estimator unchanged.}{...}
{* *! performance guide 1.0 28sep2026}{...}
{title:cnlsddf_performance -- Solver, constraints, time and memory}

{p 4 4 2}
Practical companion to {help cnlsddf_sparse}. The default remains OSQP cutting
planes with workingset(none). Performance depends on data and numerical controls.

{title:Compare settings on representative data}

{p 4 4 2}
Compare Clarabel cutbatch(1) and cutbatch(4) with workingset(none), then
once/bounded at fixed settings when restricted QPs become large. On small
samples, compare full as well. Full stores all N(N-1) ordered pairs. Cutting
planes do not guarantee lower peak memory or avoid quadratic growth in all cases.
Do not remove duplicate observations without preserving their frequency weights.

{title:Compatible options and budgets}

{p 4 4 2}
Nondefault cutbatch(), innercap(), innerpolicy() and workingset(once|bounded)
require solver(clarabel) without full. For full or OSQP comparisons, remove these
nondefault settings. fullstart(cold|interior) is for OSQP; Clarabel accepts only
the default fullstart(auto). See {help cnlsddf_sparse} for exact syntax.

{p 4 4 2}
maxiter() controls OSQP iterations per call. ipmmaxiter() controls Clarabel IPM
iterations per QP. maxcut() counts QP calls including strict refinement and
recovery. These are iteration budgets, not elapsed-second limits. The time and
memory caps in the research tables belong to external supervisors; the Stata
command does not expose those caps as options. Sampled memory supervision is
not a hard allocation limit and can overshoot before interruption.

{p 4 4 2}
Keep tolerance()/normtolerance() fixed when comparing speed. Leave adaptive
and innercap(1e-3) as the reference intermediate policy; compare monotone or a
lower innercap separately on representative data. A strict final Solved status
and the unchanged global shape/KKT certificate are always required. Never
count AlmostSolved or an uncertified finite iterate as success.

{title:Measure comparable quantities}

{p 4 4 2}
Use the same complete-case sample, row order, directions and numerical controls.
Run sequentially under a recorded thread environment, repeat and alternate the
configuration order. Same-data timing repetitions are not independent datasets.
Inspect command return codes first, then e(converged) and e(stop_reason). After
an error, previous e() values may remain and must not be treated as new results.

{p 4 4 2}
A Stata timer around cnlsddf_sparse includes Python startup and transfer;
e(seconds) measures the Python fit. Timeouts give an observed stopping time,
not a completed-fit time. Compute speed ratios only for matched, jointly
certified runs with the same budgets. Report excluded pairs and success counts.

{p 4 4 2}
e(rss) is the residual sum of squares, not memory. Operating-system peak RSS
is maximum resident memory. Python peak RSS does not include Stata memory.
matrix_payload_bytes and ipm_matrix_bytes count arrays, not process peak memory
or all factorization storage. A terminated run's last sampled RSS is a lower
bound. Stata timers do not measure memory.

{p 4 4 2}
For directional interpretation see {help cnlsddf_applied}.
