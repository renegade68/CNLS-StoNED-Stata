{smcl}
{* Public packaging documentation revision 25oct2026; estimator unchanged.}{...}
{* Interface documentation review 29sep2026; inherited by display-only dev11 candidate.}{...}
{* *! version 0.5.3-dev11  29sep2026 (skewness display correction)}{...}
{* Context estimates are point estimates; covariance diagnostics are not validated inference.}{...}
{viewerjumpto "Syntax" "convnls##syntax"}{...}
{viewerjumpto "Description" "convnls##description"}{...}
{viewerjumpto "Options" "convnls##options"}{...}
{viewerjumpto "Remarks" "convnls##remarks"}{...}
{viewerjumpto "Examples" "convnls##examples"}{...}
{viewerjumpto "Stored results" "convnls##results"}{...}
{viewerjumpto "References" "convnls##references"}{...}
{viewerjumpto "Author" "convnls##author"}{...}

{title:Title}

{phang}
{bf:convnls} {hline 2} Convex Nonparametric Least Squares (CNLS)

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:convnls} {it:depvar} {ifin}{cmd:,}
{cmdab:in:puts(}{it:varlist}{cmd:)}
[{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{cmdab:in:puts(}{it:varlist}{cmd:)}}input variables{p_end}
{syntab:Model}
{synopt:{cmdab:front:ier(}{it:string}{cmd:)}}frontier type: {bf:production} (default) or {bf:cost}{p_end}
{synopt:{opt context(varlist)}}contextual variables with common one-stage CNLS coefficients{p_end}
{synopt:{opt vce(string)}}legacy context covariance diagnostic: {bf:none} (default), {bf:ols}, or {bf:robust}; not validated inference{p_end}
{synopt:{opt rts(string)}}returns to scale: {bf:vrs} (default) or {bf:crs}{p_end}
{synopt:{opt nomon:otone}}suppress monotonicity constraints{p_end}
{syntab:Output variables}
{synopt:{cmdab:gen:erate(}{it:newvar}{cmd:)}}name for residual variable (default: {it:_convnls_resid}){p_end}
{synopt:{cmdab:fit:ted(}{it:newvar}{cmd:)}}name for fitted-value variable (default: {it:_convnls_fitted}){p_end}
{syntab:Computational}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}convergence tolerance (default: 1e-4){p_end}
{synopt:{opt rhopolicy(string)}}ADMM penalty balancing: {bf:normalized} (default) or dev4 {bf:raw}{p_end}
{synopt:{opt rho(#)}}initial ADMM penalty parameter; default -1 lets {cmd:convnls} choose automatically{p_end}
{synopt:{cmdab:max:iter(}{it:#}{cmd:)}}maximum ADMM iterations per inner call (default: 20000){p_end}
{synopt:{cmdab:maxc:ut(}{it:#}{cmd:)}}maximum cutting-plane outer iterations (default: 100){p_end}
{synopt:{opt addper(#)}}maximum number of newly added violations per cutting-plane iteration under {cmd:vcstrat(all)}; 0 means no user cap (internal candidate batch limit: 200000){p_end}
{synopt:{opt prune(#)}}prune constraints with slack >= {it:#} (default: 0, no pruning){p_end}
{synopt:{opt level(#)}}accepted for compatibility; default is Stata's current confidence level; no confidence intervals are produced{p_end}
{synopt:{opt verbose}}display iteration-level progress{p_end}
{syntab:Algorithm (CNLS+-G)}
{synopt:{opt init(string)}}constraint initialisation: {bf:none} (default), {bf:ele}, {bf:sweet}, or {bf:adam}{p_end}
{synopt:{opt sweetpct(#)}}sweet-spot distance percentile for {bf:init(sweet)} (default: 3){p_end}
{synopt:{opt vcstrat(string)}}violation strategy: {bf:all} (default), {bf:group}, or {bf:one}{p_end}
{synopt:{opt adamiter(#)}}Adam warm-start iterations for {bf:init(adam)} (default: 100){p_end}
{synopt:{opt adamlr(#)}}Adam warm-start learning rate (default: 0.01){p_end}
{synopt:{opt adamlambda(#)}}penalty weight used in the Adam warm start (default: 10){p_end}
{synopt:{opt adamscreen(#)}}maximum number of Adam-screened candidate constraints; 0 means no user cap (internal candidate batch limit: 200000) (default: 2000){p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:convnls} fits the Convex Nonparametric Least Squares (CNLS) estimator of
Kuosmanen (2008).  CNLS is a shape-constrained regression that imposes
piecewise-linear convexity (cost frontier) or concavity (production frontier)
on the conditional mean function, without assuming a parametric functional form.
It is the first stage of the Stochastic Non-smooth Envelopment of Data
(StoNED) estimator (Kuosmanen and Kortelainen, 2012).

{pstd}
The estimator solves the quadratic programme

{pmore}
min  sum_i (y_i - alpha_i - beta_i'x_i)^2

{pmore}
subject to: alpha_i + beta_i'x_i >= alpha_h + beta_h'x_i  for all i != h
(Afriat inequalities, cost frontier; reversed for production)

{pmore}
and  beta_i >= 0  for all i  (monotonicity, if not suppressed)

{pstd}
where {it:alpha_i} is an observation-specific intercept and
{it:beta_i} is an observation-specific gradient vector.
Under variable returns to scale ({opt rts(vrs)}) these parameters are free.
Under constant returns to scale ({opt rts(crs)}) alpha_i = 0 for all i.
In additive StoNED, CRS needs some care: as discussed by Kuosmanen and
Kortelainen (2012, Section 4.2), imposing CRS on the first-stage
average-practice CNLS function and then shifting by expected inefficiency
does not generally guarantee that the shifted best-practice frontier also
satisfies CRS.  Multiplicative or log specifications are the natural
setting when scale homogeneity is central, but need their own model restrictions.
The command applies convexity/concavity and {opt rts(crs)} to the variables
actually supplied. In particular, alpha_i=0 in log coordinates does not impose
degree-one homogeneity in original quantities. Supplying logged variables is
not the separate multiplicative CNLS optimisation problem.

{pstd}
The CNLS residuals epsilon_i = y_i - (alpha_i + beta_i'x_i) are saved and
passed to the second stage.  In the StoNED framework, the second stage
decomposes the residuals into a symmetric noise component (v) and a
one-sided inefficiency component (u) using {help stoned}.

{pstd}
Residuals and fitted values are returned in the supplied dependent variable's
units. The commands do not automatically log-transform outcomes. The exponential
score from {cmd:stoned} is a conditional efficiency ratio under a suitable
natural-log error model; for a dependent variable in levels it is a unit-dependent
index. The conditional mean of u remains in residual units. See
{help stoned##units:Residual units and interpretation}.

{pstd}
With {opt context()}, {cmd:convnls} estimates common linear coefficients for
contextual variables jointly with an additive CNLS frontier in the coordinates
supplied by the user. This is related to the one-stage approach of Johnson and
Kuosmanen (2011), but does not automatically implement their multiplicative
level-frontier specification by supplying a logged dependent variable.
Only point estimates are supported. This version does not post {cmd:e(V)} or
report standard errors, p-values, or confidence intervals for context coefficients,
the nonparametric frontier, or subsequent StoNED efficiency scores.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{cmdab:in:puts(}{it:varlist}{cmd:)} specifies the input (regressor) variables.
At least one input variable is required.

{phang}
{cmdab:front:ier(}{it:string}{cmd:)} specifies the frontier orientation.
{opt frontier(production)} (the default) imposes concavity.
{opt frontier(cost)} imposes convexity. First-stage CNLS estimates a
conditional-mean function; its fitted values need not envelop all observations.

{phang}
{opt context(varlist)} specifies contextual variables whose common linear
coefficients are estimated jointly with the CNLS frontier.  The option is
currently supported with {opt rts(vrs)}. The coefficients are stored in {cmd:e(b)}.
The dependent variable cannot be included in {cmd:context()}.
Constant or numerically dependent context columns and affine aliases with the
inputs are rejected before writing outputs. The check centers and normalizes
columns and uses a relative singular-value threshold of 1e-10. Passing this
necessary design check does not establish general identification or statistical
inference. A nonconverged iterate is provisional; inspect {cmd:e(converged)}.

{phang}
{opt vce(string)} now selects an optional legacy covariance diagnostic, not an
inferential VCE. The default is {opt vce(none)}. Explicit {opt vce(ols)} or
{opt vce(robust)} saves {cmd:e(V_context_diagnostic)} and
{cmd:e(V_robust_diagnostic)} only after strict or practical CNLS certification.
With no {opt context()}, these options have no effect and post no diagnostics.
The old names {cmd:e(V)} and {cmd:e(V_robust)} are not posted. This is an
intentional compatibility change from dev9 and the accepted implementation.

{pstd}
For centered context matrix C, H=(C'C)^(-1), and CNLS residual vector r, the
legacy OLS and sandwich formulas are RSS/(N-q)*H and
N/(N-q)*H*C'diag(r^2)*C*H. {cmd:e(V_context_diagnostic)} holds the formula
selected by {opt vce()}; {cmd:e(V_robust_diagnostic)} always holds the sandwich
formula. Under {opt vce(robust)}, the two matrices are identical.
They are retained for numerical comparison, not as established asymptotic
covariances of this implementation. N-q is not an established effective residual
degrees of freedom for the jointly fitted CNLS frontier. A sandwich calculation
does not by itself correct frontier bias, nonidentification, or dependence.
The {opt level()} option remains accepted for syntax compatibility but produces
no confidence intervals. No claim of interval coverage is made.

{phang}
{opt rts(string)} specifies the returns-to-scale assumption.
{opt rts(vrs)} (the default) allows variable returns to scale by including
free firm-specific intercepts {it:alpha_i}.
{opt rts(crs)} imposes constant returns to scale by constraining {it:alpha_i = 0}.
This restriction is applied to the first-stage average-practice function;
interpret the shifted StoNED frontier under CRS in light of the additive
versus multiplicative distinction noted above.

{phang}
{opt nomonotone} suppresses the input monotonicity constraints
{it:beta_ij >= 0}.  By default, monotonicity is imposed: production (cost)
output is non-decreasing (non-decreasing) in each input.

{dlgtab:Output variables}

{phang}
{cmdab:gen:erate(}{it:newvar}{cmd:)} specifies the name of the new variable
to hold the CNLS residuals.  If not specified, the residuals are stored in
{it:_convnls_resid}.

{phang}
{cmdab:fit:ted(}{it:newvar}{cmd:)} specifies the name of the new variable to
hold the CNLS fitted values.  If not specified, the fitted values are stored
in {it:_convnls_fitted}.

{dlgtab:Computational}

{phang}
{cmdab:tol:erance(}{it:#}{cmd:)} sets the ADMM primal/dual convergence
tolerance.  Default is 1e-4.  Tighten to 1e-5 or 1e-6 for more precise
residuals when n is small.  A run is reported as practically converged
({cmd:e(converged)=2}) when the maximum Afriat violation remains below
30 times this tolerance after two consecutive checks with certified inner
optimality and monotonicity. No success status is assigned to an inner iteration
limit or numerical failure.

{phang}
{opt rho(#)} sets the initial ADMM penalty parameter.  The default -1 lets
{cmd:convnls} choose the starting value automatically.  This is an advanced
diagnostic option and is not usually needed in routine use.

{phang}
{cmdab:max:iter(}{it:#}{cmd:)} sets the maximum number of ADMM inner
iterations per inner solver call. Default is 20000. A cutting-plane step can
include an additional tightening call, with its own budget. Calls may continue
saved ADMM state when the working set is unchanged; {opt maxiter()} is not a
limit on total ADMM iterations.

{phang}
{cmdab:maxc:ut(}{it:#}{cmd:)} sets the maximum number of cutting-plane
outer iterations. Default is 100. With {opt vcstrat(one)}, a value of 100
is automatically replaced by min(N*(N-1),5000), including an explicitly supplied
{opt maxcut(100)}. Each step solves or continues the restricted ADMM problem
and scans all Afriat pairs. Finding no new violations is not sufficient to
terminate successfully; inner optimality and global shape certification must
also pass.

{phang}
{opt prune(#)} removes constraints whose slack is at least {it:#}, in internal
scaled coordinates, before solves starting at outer iteration 4.  A value such as {opt prune(1e-4)} can speed
up computations for large samples by reducing the active constraint set.
Default is 0 (no pruning).
Slack-only pruning can repeatedly remove and re-add constraints and exhaust
{opt maxcut()}. Check the convergence code; use {opt prune(0)} when diagnosing
such behavior. A smaller RSS from an uncertified pruned fit is not an improvement.

{phang}
{opt verbose} requests iteration-level output showing the maximum
constraint violation and ADMM residuals at each cutting-plane step.

{dlgtab:Algorithm (CNLS+-G)}

{phang}
{opt init(string)} selects the initial set of Afriat constraints before the
cutting-plane loop begins.
{p_end}
{pmore}
{opt init(none)} (the default) starts the cutting-plane loop with no Afriat
constraints, which is the standard approach.
{p_end}
{pmore}
{opt init(ele)} applies the Elementary Afriat initialisation of
Lee et al. (2013): observations are sorted by the first input variable and
the 2(n-1) adjacent-pair constraints are added as the initial set.
This provides a structured warm start that reduces the number of
cutting-plane iterations in most samples.
{p_end}
{pmore}
{opt init(sweet)} extends {opt init(ele)} with the Sweet Spot approach of
Lee et al. (2013): for each observation i, all pairs (i,h) for which the
Euclidean distance d(i,h) in internally scaled input space is within the p-th percentile of
distances from i are added with their reverse pairs (zero cutoffs add no
extra neighbors), where p is controlled by
{opt sweetpct()}.  These constraints are solved in the first cutting-plane iteration before
any convergence check.
This is the most aggressive initialisation and is recommended for n >= 100.
{p_end}
{pmore}
{opt init(adam)} uses a short penalized Adam warm start to screen likely
violated Afriat constraints before the deterministic cutting-plane solve.
The Adam step only seeds the active set; the reported CNLS estimates are
still obtained by solving the constrained least-squares problem with final
violation checks.
{p_end}

{phang}
{opt sweetpct(#)} sets the distance percentile threshold for the Sweet Spot
initialisation.  Default is 3 (3rd percentile).  Higher values include more
initial constraints but may increase the size of the initial solve.
Only relevant when {opt init(sweet)} is specified.
{p_end}

{phang}
{opt vcstrat(string)} selects the strategy for choosing which violated Afriat
constraints to add at each cutting-plane iteration.
{p_end}
{pmore}
{opt vcstrat(all)} (the default) selects newly violated constraints from a complete pairwise scan.
With {opt addper(0)}, candidate storage retains at most the 200000 largest
new violations per iteration; a positive {opt addper()} sets the candidate cap
instead (up to N*(N-1)).
the scan always continues to compute the true global maximum violation.
Use {opt addper(#)} to cap the number of newly added constraints to the
largest {it:#} violations per iteration.
{p_end}
{pmore}
{opt vcstrat(group)} implements the CNLS-G strategy of Lee et al. (2013):
for each observation i, only the most violated currently inactive constraint for that i
is added per iteration. Already active violations remain in the global diagnostics.  At most n constraints are added per step.
Each step is faster but more iterations are needed.
This is recommended in combination with {opt init(sweet)} for large samples.
{p_end}
{pmore}
{opt vcstrat(one)} implements the CNLS-O strategy: only the single globally
most violated inactive (i,h) pair is added per iteration.  This uses the fewest
constraints per step but requires the most iterations overall, and is
typically dominated by {opt vcstrat(group)}.
{p_end}

{phang}
{ul:Convergence certification.}
All combinations of {opt init()} and {opt vcstrat()} target the same CNLS
least-squares problem. The restricted QP is solved before every convergence
check, including the initial check. A success code requires certified inner
primal, dual, stationarity and complementarity diagnostics, monotonicity, and
an exhaustive Afriat check. Different strategies can still return different
approximations at finite tolerance or when an iteration budget is exhausted.
Inspect {cmd:e(converged)} and {cmd:e(stop_reason)}; a feasible iterate alone is
not evidence of an optimum.
{p_end}


{marker remarks}{...}
{title:Remarks}

{pstd}
{ul:Algorithm.}
{cmd:convnls} uses a cutting-plane outer loop combined with an
Alternating Direction Method of Multipliers (ADMM) inner solver.
By default, the cutting-plane loop starts with no Afriat constraints and
iteratively adds newly violated constraints from a complete pairwise scan.
Termination also requires a certified restricted QP solve.
Constraint deduplication prevents duplicate pairs in the current working set.
A pruned pair may be added again.
Optional CNLS+-G acceleration (Lee et al. 2013) is available through the
{opt init()} and {opt vcstrat()} options: {opt init(sweet)} pre-populates
the active set with nearby-observation pairs before the loop begins, and
{opt vcstrat(group)} adds only the worst violation per observation per step.
These options do not change the target CNLS problem; runtime gains depend on
the data and settings.

{pstd}
{ul:Scaling.}
The dependent variable is centered and scaled under VRS, and only scaled
under CRS. Inputs are divided by their standard deviations without centering
to improve numerical conditioning.  Contextual variables specified in {opt context()} are centered
and standardised internally, and their reported coefficients are transformed
back to the original scale.  Fitted values and residuals are returned on the
original scale.

{pstd}
{ul:Saved variables.}
In addition to the variables named by {opt generate()} and {opt fitted()},
{cmd:convnls} saves {it:_convnls_alpha} (intercepts) and
{it:_convnls_beta1}, {it:_convnls_beta2}, ... (gradients) in the dataset.
These are overwritten if {cmd:convnls} is called again.
Output names must be distinct and must not coincide with the dependent variable,
inputs, contextual variables or coefficient destinations. Matching uses exact
names, not variable abbreviations. Existing exact-name output variables are
replaced on normal completion, including an uncertified result explicitly marked
by {cmd:e(converged)=0}. Invalid inputs or a runtime error leave the dataset as it
was before the call; an error is not a new successful estimate. Data protection
uses {cmd:preserve} and may require additional memory or temporary disk space.

{pstd}
{ul:Residual skewness display.}
Undefined skewness (including a constant residual) is reported without a sign
comparison. Exactly zero skewness is reported as having no direction. Finite
nonzero skewness keeps its sign, however small; no near-zero threshold is applied.
A sign-consistency note is descriptive and does not certify a valid variance
decomposition, convergence, or distributional assumptions.

{pstd}
{ul:StoNED workflow.}
The standard two-stage StoNED estimator uses {help convnls} (Stage 1)
followed by {help stoned} (Stage 2):

{phang2}(1) Estimate CNLS and save residuals:{p_end}
{phang3}{cmd:convnls y, inputs(x1 x2) frontier(cost) generate(ehat)}{p_end}

{phang2}(2) Decompose residuals via method of moments and compute efficiency scores:{p_end}
{phang3}{cmd:stoned ehat, implementation(legacy) frontier(cost) generate(te)}{p_end}

{pstd}
{cmd:stoned} estimates sigma_u and sigma_v by the method of moments
(Kuosmanen and Kortelainen, 2012) and computes firm-level efficiency
scores using the Battese and Coelli (1988) conditional expectation formula.
Using Stata's built-in {cmd:frontier} on CNLS residuals is not recommended:
first-stage CNLS residuals need not retain the mean structure of a
normal--half-normal composed error, causing the MLE to converge to the
boundary sigma_u=0.

{pstd}
{ul:Computational cost.}
The number of potential Afriat constraints is N(N-1), which grows as O(N^2).
Writing P=N*(1+m)+q under VRS and A for the current number of constraints,
the structured engine stores Afriat observation pairs and their local coefficient
blocks instead of dense A-by-P constraint matrices. Constraint products cost
O(A*(1+m)+P), and the Gram matrix is assembled from local blocks. Cholesky
factors remain cached within each restricted QP solve, with O(P^3)
factorization after a working-set or penalty change and O(P^2) triangular
solves per iteration. Normal, Gram, system and factor matrices remain dense
P-by-P matrices. The complete Afriat scan still costs O(N^2*m).
{cmd:e(constraint_bytes)} reports the largest compact operator matrix payload;
{cmd:e(dense_constraint_bytes)} reports the largest equivalent single dense
constraint matrix payload. These exclude temporary arrays, allocator overhead,
other model matrices and the observation-pair tracking matrix; neither measures
peak process memory. Summation order can change numerical trajectories, but
objective, scaling and convergence thresholds are unchanged.
For larger samples, estimation may take several minutes.
Use {opt verbose} to monitor progress.

{pstd}
{ul:Convergence status.}
{cmd:e(converged)} is 0 for an uncertified result, 1 when the maximum Afriat
and monotonicity violations are at most {opt tolerance()}, and 2 when the
Afriat violation is below 30*{opt tolerance()} in two successive certified
checks. Both success codes also require inner optimality diagnostics at an
inner tolerance no larger than the requested tolerance and monotonicity at
the requested tolerance. Shape violations are reported in the internally
scaled coordinates, including positive violations below the tolerance.
Diagnostics refer to the final returned iterate, including budget exits.
{cmd:e(stop_reason)} distinguishes convergence, practical convergence,
outer iteration limits, inner iteration limits, and numerical failure.


{pstd}
When the ordered working constraints are unchanged, repeated inner solves retain
ADMM slack and scaled dual variables, penalty, relaxation and iteration schedule.
This also applies when tightening the inner tolerance. Adding or pruning a
constraint discards this state; a new estimation always starts fresh. All
convergence checks are recomputed at the current tolerance. The matrix factor
cache is local to each inner call and is rebuilt on continuation.

{pstd}
With {cmd:rhopolicy(normalized)}, the adaptive penalty compares each primal/dual residual divided by its existing
stopping threshold. Every 100 cumulative ADMM iterations before iteration 5,000
of the same restricted QP, a ratio imbalance above
5 doubles or halves the penalty, with adaptive bounds 1e-6 and 1e6. The scaled dual
variable is rescaled to preserve the unscaled multiplier. Initial {cmd:rho()} values
outside those bounds are accepted; an update never reverses direction to reach a
bound. From iteration 5,000 onward the penalty is held fixed to avoid sustained
oscillation from repeated penalty changes. ADMM and all certification checks
continue within the requested budgets. Continuation (including a tighter tolerance)
preserves this clock; rebuilding the working set resets it. These ratios guide
penalty updates only. Objective, data scaling, primal,
dual, stationarity, complementarity and global shape certification are unchanged.
This policy is a heuristic and does not guarantee faster convergence for every sample.
Use {cmd:rhopolicy(raw)} to reproduce dev4 penalty balancing: raw residuals,
its original 0.1 floor, no added upper bound and no adaptation cutoff.
This can be faster on some samples.
The reported {cmd:e(converged)} still distinguishes strict (1), practical (2)
and uncertified (0) termination under either policy.

{marker examples}{...}
{title:Examples}

{pstd}
The simulated examples below are self-contained. Save any data in memory
before running them. A short scalar example is included in examples/quickstart.do.
Empirical bank data and manuscript replication scripts are distributed separately.

{pstd}{ul:Example 1: Simulated cost frontier (VRS), StoNED two-stage}

{phang2}{stata clear:. clear}{p_end}
{phang2}{stata set seed 42:. set seed 42}{p_end}
{phang2}{stata set obs 40:. set obs 40}{p_end}
{phang2}{stata gen ln_x1 = rnormal(1.5, 0.6):. gen ln_x1 = rnormal(1.5, 0.6)}{p_end}
{phang2}{stata gen ln_x2 = rnormal(1.0, 0.4):. gen ln_x2 = rnormal(1.0, 0.4)}{p_end}
{phang2}{stata gen u = abs(rnormal(0, 0.20)):. gen u = abs(rnormal(0, 0.20))}{p_end}
{phang2}{stata gen ln_c = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) + u:. gen ln_c = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) + u}{p_end}

{phang2}{stata convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat) fitted(fhat):. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat) fitted(fhat)}{p_end}
{phang2}{stata stoned ehat, implementation(legacy) frontier(cost) generate(te) inefficiency(uhat):. stoned ehat, implementation(legacy) frontier(cost) generate(te) inefficiency(uhat)}{p_end}
{phang2}{stata summarize te uhat:. summarize te uhat}{p_end}

{pstd}{ul:Example 2: Simulated production frontier (CRS)}

{phang2}{stata clear:. clear}{p_end}
{phang2}{stata set seed 4242:. set seed 4242}{p_end}
{phang2}{stata set obs 100:. set obs 100}{p_end}
{phang2}{stata gen ln_x1 = rnormal(1.5, 0.6):. gen ln_x1 = rnormal(1.5, 0.6)}{p_end}
{phang2}{stata gen ln_x2 = rnormal(1.0, 0.4):. gen ln_x2 = rnormal(1.0, 0.4)}{p_end}
{phang2}{stata gen u = abs(rnormal(0, 0.20)):. gen u = abs(rnormal(0, 0.20))}{p_end}
{phang2}{stata gen ln_y = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) - u:. gen ln_y = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) - u}{p_end}

{phang2}{stata convnls ln_y, inputs(ln_x1 ln_x2) frontier(production) rts(crs) generate(ehat_p) init(sweet) sweetpct(3) vcstrat(group) maxcut(500):. convnls ln_y, inputs(ln_x1 ln_x2) frontier(production) rts(crs) generate(ehat_p) init(sweet) sweetpct(3) vcstrat(group) maxcut(500)}{p_end}
{phang2}{stata stoned ehat_p, implementation(legacy) frontier(production) generate(te_p):. stoned ehat_p, implementation(legacy) frontier(production) generate(te_p)}{p_end}

{pstd}{ul:Example 3: Constant returns to scale}

{phang2}{stata clear:. clear}{p_end}
{phang2}{stata set seed 42:. set seed 42}{p_end}
{phang2}{stata set obs 40:. set obs 40}{p_end}
{phang2}{stata gen ln_x1 = rnormal(1.5, 0.6):. gen ln_x1 = rnormal(1.5, 0.6)}{p_end}
{phang2}{stata gen ln_x2 = rnormal(1.0, 0.4):. gen ln_x2 = rnormal(1.0, 0.4)}{p_end}
{phang2}{stata gen u = abs(rnormal(0, 0.20)):. gen u = abs(rnormal(0, 0.20))}{p_end}
{phang2}{stata gen ln_c = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) + u:. gen ln_c = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) + u}{p_end}
{phang2}{stata convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) rts(crs) generate(ehat_crs):. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) rts(crs) generate(ehat_crs)}{p_end}

{pstd}{ul:Example 4: Verbose output and convergence diagnostics}

{phang2}{stata convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat_diag) verbose tolerance(1e-5):. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat_diag) verbose tolerance(1e-5)}{p_end}
{phang2}{stata ereturn list:. ereturn list}{p_end}

{pstd}{ul:Example 5: CNLS+-G accelerated algorithm (Lee et al. 2013)}

{phang2}{cmd:. * Recommended for n >= 100: Sweet Spot init + Group violation strategy}{p_end}
{phang2}{stata convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat_sweet) init(sweet) sweetpct(3) vcstrat(group):. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat_sweet) init(sweet) sweetpct(3) vcstrat(group)}{p_end}

{phang2}{cmd:. * Check which algorithm was used:}{p_end}
{phang2}{stata ereturn list:. ereturn list}{p_end}

{pstd}{ul:Example 6: Elementary Afriat init only}

{phang2}{stata convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat_ele) init(ele) vcstrat(group):. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat_ele) init(ele) vcstrat(group)}{p_end}

{pstd}{ul:Example 7: Joint contextual-variable point estimates}

{phang2}{stata clear:. clear}{p_end}
{phang2}{stata set seed 270627:. set seed 270627}{p_end}
{phang2}{stata set obs 60:. set obs 60}{p_end}
{phang2}{stata gen x = runiform()*10 + 1:. gen x = runiform()*10 + 1}{p_end}
{phang2}{stata gen z = rnormal():. gen z = rnormal()}{p_end}
{phang2}{stata gen u = abs(rnormal(0,.05)):. gen u = abs(rnormal(0,.05))}{p_end}
{phang2}{stata gen y = ln(x) + 2 + .35*z + rnormal(0,.05) - u:. gen y = ln(x) + 2 + .35*z + rnormal(0,.05) - u}{p_end}

{phang2}{stata convnls y, inputs(x) context(z) frontier(production) init(sweet) vcstrat(group) vce(robust) generate(ehat_ctx) fitted(fhat_ctx):. convnls y, inputs(x) context(z) frontier(production) init(sweet) vcstrat(group) vce(robust) generate(ehat_ctx) fitted(fhat_ctx)}{p_end}
{phang2}{stata ereturn list:. ereturn list}{p_end}


{pstd}
{ul:Residual provenance for StoNED.}
The generated residual stores CNLS version, frontier, RTS, original sample size,
convergence and inner-solver status, stop reason, and a value signature as
variable characteristics ({cmd:cnls_*}). These survive {cmd:save}/{cmd:use},
{cmd:rename} and a full {cmd:clonevar}. An uncertified final iterate is still
returned by {cmd:convnls}, with convergence 0 recorded on the residual;
{help stoned} refuses its use. Strict (1) and practical (2) certification remain
distinct and both are supported by {cmd:stoned}.

{pstd}
The signature checks the residual-value multiset, including missing values and
observation count. It is independent of variable name, row order and display
format, and is computed on a double copy without sorting the dataset. It does
not authenticate metadata, link each residual to other variables, or check
whether the original inputs/outcome have since changed. Copying values with
{cmd:generate} or through CSV does not carry these characteristics.

{pstd}
{ul:Estimation sample.}
On normal completion, {cmd:e(sample)} identifies the observations actually used,
with or without {opt context()}; their count equals {cmd:e(N)}. The sample obeys
{cmd:if}/{cmd:in} and excludes missing dependent-variable values (including
extended missing values). Missing inputs or context variables within that candidate
sample cause r(459), preserving the existing policy; they are not silently omitted.
Missing predictors outside the candidate sample do not cause this error.
Generated residuals, fitted values and observation-specific coefficients are
missing outside the estimation sample.

{pstd}
For an uncertified final iterate, {cmd:e(sample)} still describes the attempted
fit and is not a convergence certificate. Check {cmd:e(converged)} as well.
A later estimation replaces the current sample; before {cmd:stoned}, retain the
first-stage mask with {cmd:generate byte cnls_sample = e(sample)} if needed.
Without context, posting the sample does not create a parametric coefficient
vector or covariance matrix. No new predict or inference support is implied.

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:convnls} stores the following in {cmd:e()}, including the function {cmd:e(sample)}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(m)}}number of input variables{p_end}
{synopt:{cmd:e(q)}}number of contextual variables specified in {cmd:context()}{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}
{synopt:{cmd:e(sigma2)}}estimated variance (RSS / N){p_end}
{synopt:{cmd:e(converged)}}0 not converged, 1 converged, 2 practically converged{p_end}
{synopt:{cmd:e(cp_iters)}}number of cutting-plane iterations{p_end}
{synopt:{cmd:e(n_active)}}number of Afriat constraints in the final solved working set{p_end}
{synopt:{cmd:e(max_viol)}}true maximum positive Afriat violation of the returned fit, on the internal scale{p_end}
{synopt:{cmd:e(max_mono)}}maximum positive monotonicity violation, internal scale; 0 when disabled{p_end}
{synopt:{cmd:e(inner_status)}}last inner solve: 1 certified, 0 iteration limit, -1 numerical failure (including an uncertified unconstrained solve){p_end}
{synopt:{cmd:e(inner_iter)}}iterations in the last inner solve{p_end}
{synopt:{cmd:e(inner_primal)}, {cmd:e(inner_dual)}}last inner primal and dual residual infinity norms{p_end}
{synopt:{cmd:e(inner_stationarity)}}unregularized QP stationarity infinity norm{p_end}
{synopt:{cmd:e(inner_complementarity)}}maximum absolute multiplier-slack product{p_end}
{synopt:{cmd:e(inner_eps_primal)}, {cmd:e(inner_eps_dual)}}inner residual thresholds{p_end}
{synopt:{cmd:e(inner_eps_stationarity)}, {cmd:e(inner_eps_complementarity)}}inner optimality thresholds{p_end}
{synopt:{cmd:e(inner_tolerance)}}tolerance used by the last inner solve{p_end}
{synopt:{cmd:e(total_inner_iters)}}total ADMM iterations across all solves{p_end}
{synopt:{cmd:e(inner_solves)}}number of restricted QP solves, including tightening solves{p_end}
{synopt:{cmd:e(warm_starts)}}inner calls continuing saved ADMM state for unchanged constraints{p_end}
{synopt:{cmd:e(rho_adapt_limit)}}iteration at which penalty adaptation stops: 5000 for normalized; 0 means unlimited for raw{p_end}
{synopt:{cmd:e(final_rho)}}ADMM penalty used for the returned iterate{p_end}
{synopt:{cmd:e(factorizations)}}total Cholesky factorization attempts across inner solves{p_end}
{synopt:{cmd:e(chol_solves)}}total cached-factor triangular solve attempts{p_end}
{synopt:{cmd:e(lu_fallbacks)}}total LU fallback attempts after Cholesky failure{p_end}
{synopt:{cmd:e(rho_updates)}}total ADMM penalty changes requiring a new factorization{p_end}
{synopt:{cmd:e(constraint_bytes)}}largest structured constraint-operator matrix payload in bytes{p_end}
{synopt:{cmd:e(dense_constraint_bytes)}}largest equivalent single dense constraint-matrix payload in bytes{p_end}
{synopt:{cmd:e(addper)}}cap on newly added violations per iteration; 0 means no user cap (internal candidate batch limit: 200000){p_end}
{synopt:{cmd:e(adamiter)}}Adam warm-start iterations{p_end}
{synopt:{cmd:e(adamlr)}}Adam warm-start learning rate{p_end}
{synopt:{cmd:e(adamlambda)}}Adam warm-start penalty weight{p_end}
{synopt:{cmd:e(adamscreen)}}Adam candidate-screening cap{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(rho_policy)}}penalty balancing rule: normalized or raw{p_end}
{synopt:{cmd:e(stop_reason)}}converged, practical, maxcut, inner_maxiter, or inner_failure{p_end}
{synopt:{cmd:e(engine)}}{cmd:structured}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:convnls}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(inputs)}}names of input variables{p_end}
{synopt:{cmd:e(context)}}names of contextual variables, if specified{p_end}
{synopt:{cmd:e(frontier)}}{cmd:production} or {cmd:cost}{p_end}
{synopt:{cmd:e(rts)}}{cmd:vrs} or {cmd:crs}{p_end}
{synopt:{cmd:e(version)}}CNLS implementation version{p_end}
{synopt:{cmd:e(residual_schema)}}residual provenance schema ({cmd:1}){p_end}
{synopt:{cmd:e(residual_signature)}}signature recorded on the generated residual{p_end}
{synopt:{cmd:e(resid_var)}}name of the residual variable created{p_end}
{synopt:{cmd:e(fitted_var)}}name of the fitted-value variable created{p_end}
{synopt:{cmd:e(init_method)}}{cmd:none}, {cmd:ele}, {cmd:sweet}, or {cmd:adam}; initialisation used{p_end}
{synopt:{cmd:e(vcstrat)}}{cmd:all}, {cmd:group}, or {cmd:one}; violation strategy used{p_end}
{synopt:{cmd:e(vce)}}{cmd:none}; no inferential covariance is posted{p_end}
{synopt:{cmd:e(inference_status)}}{cmd:not_validated} with {cmd:context()}{p_end}
{synopt:{cmd:e(context_diagnostic)}}requested legacy diagnostic with {cmd:context()}: none, ols, or robust{p_end}
{synopt:{cmd:e(context_diagnostic_status)}}with {cmd:context()}: not_requested, not_certified, or legacy_formula_only{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}contextual-variable coefficient estimates, if {cmd:context()} is specified{p_end}
{synopt:{cmd:e(V_context_diagnostic)}}selected legacy formula, only when explicitly requested and certified; not {cmd:e(V)}{p_end}
{synopt:{cmd:e(V_robust_diagnostic)}}legacy sandwich formula under the same restrictions; not validated inference{p_end}
{synopt:{cmd:e(cp_history)}}cutting-plane iteration history: iteration, solved working-set size, total violations above tolerance (including active pairs), newly added constraints, true global maximum violation{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}
Badunenko, O., P. Mozharovskyi, and Y. Kolomiytseva. 2020.
{it:npsf: Nonparametric and stochastic efficiency and productivity analysis}.
R package version 0.8.0.
{browse "https://CRAN.R-project.org/package=npsf":CRAN package page}.
DOI: {browse "https://doi.org/10.32614/CRAN.package.npsf":10.32614/CRAN.package.npsf}.{p_end}

{phang}
Badunenko, O., and P. Mozharovskyi. 2016.
Nonparametric frontier analysis using Stata.
{it:Stata Journal} 16(3): 550-589.
DOI: {browse "https://doi.org/10.1177/1536867X1601600302":10.1177/1536867X1601600302}.{p_end}

{phang}
Kuosmanen, T. 2008.
Representation theorem for convex nonparametric least squares.
{it:Econometrics Journal} 11(2): 308-325.
DOI: {browse "https://doi.org/10.1111/j.1368-423X.2008.00239.x":10.1111/j.1368-423X.2008.00239.x}.{p_end}

{phang}
Kuosmanen, T., and M. Kortelainen. 2012.
Stochastic non-smooth envelopment of data: semi-parametric
frontier estimation subject to shape constraints.
{it:Journal of Productivity Analysis} 38(1): 11-28.
DOI: {browse "https://doi.org/10.1007/s11123-010-0201-3":10.1007/s11123-010-0201-3}.{p_end}

{phang}
Lee, C.-Y., A. L. Johnson, E. Moreno-Centeno, and T. Kuosmanen. 2013.
A more efficient algorithm for convex nonparametric least squares.
{it:European Journal of Operational Research} 227(2): 391-400.
DOI: {browse "https://doi.org/10.1016/j.ejor.2012.11.054":10.1016/j.ejor.2012.11.054}.{p_end}


{marker author}{...}
{title:Author}

{pstd}
Choonjoo Lee{break}
Department of Defense AI/Robotics, Korea National Defense University{break}
Nonsan, Republic of Korea{break}
sarang90@korea.kr

{pstd}
Bug reports and suggestions are welcome.
