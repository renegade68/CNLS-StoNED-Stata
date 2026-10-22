{smcl}
{* Help revision: 2026-10-22.1; computational version unchanged.}{...}
{* *! version 0.5.2  07jul2026 (help examples and references verified)}{...}
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
{synopt:{cmdab:fr:ontier(}{it:string}{cmd:)}}frontier type: {bf:production} (default) or {bf:cost}{p_end}
{synopt:{opt context(varlist)}}contextual variables with common one-stage CNLS coefficients{p_end}
{synopt:{opt vce(string)}}variance estimator for {cmd:context()} coefficients: {bf:ols} (default) or {bf:robust}{p_end}
{synopt:{opt rts(string)}}returns to scale: {bf:vrs} (default) or {bf:crs}{p_end}
{synopt:{opt nomo:notone}}suppress monotonicity constraints{p_end}
{syntab:Output variables}
{synopt:{cmdab:gen:erate(}{it:newvar}{cmd:)}}name for residual variable (default: {it:_convnls_resid}){p_end}
{synopt:{cmdab:fit:ted(}{it:newvar}{cmd:)}}name for fitted-value variable (default: {it:_convnls_fitted}){p_end}
{syntab:Computational}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}convergence tolerance (default: 1e-4){p_end}
{synopt:{opt rho(#)}}initial ADMM penalty parameter; default -1 lets {cmd:convnls} choose automatically{p_end}
{synopt:{cmdab:maxi:ter(}{it:#}{cmd:)}}maximum ADMM iterations (default: 20000){p_end}
{synopt:{cmdab:maxc:ut(}{it:#}{cmd:)}}maximum cutting-plane outer iterations (default: 100){p_end}
{synopt:{opt addper(#)}}maximum number of newly added violations per cutting-plane iteration under {cmd:vcstrat(all)}; 0 means no cap{p_end}
{synopt:{opt prune(#)}}prune constraints with slack >= {it:#} (default: 0, no pruning){p_end}
{synopt:{opt verbose}}display iteration-level progress{p_end}
{syntab:Algorithm (CNLS+-G)}
{synopt:{opt init(string)}}constraint initialisation: {bf:none} (default), {bf:ele}, {bf:sweet}, or {bf:adam}{p_end}
{synopt:{opt sweetpct(#)}}sweet-spot distance percentile for {bf:init(sweet)} (default: 3){p_end}
{synopt:{opt vcstrat(string)}}violation strategy: {bf:all} (default), {bf:group}, or {bf:one}{p_end}
{synopt:{opt adamiter(#)}}Adam warm-start iterations for {bf:init(adam)} (default: 100){p_end}
{synopt:{opt adamlr(#)}}Adam warm-start learning rate (default: 0.01){p_end}
{synopt:{opt adamlambda(#)}}penalty weight used in the Adam warm start (default: 10){p_end}
{synopt:{opt adamscreen(#)}}maximum number of Adam-screened candidate constraints; 0 means no cap (default: 2000){p_end}
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
setting when scale homogeneity is central.

{pstd}
The CNLS residuals epsilon_i = y_i - (alpha_i + beta_i'x_i) are saved and
passed to the second stage.  In the StoNED framework, the second stage
decomposes the residuals into a symmetric noise component (v) and a
one-sided inefficiency component (u) using {help stoned}.

{pstd}
With {opt context()}, {cmd:convnls} estimates common linear coefficients for
contextual variables jointly with the CNLS frontier, following the one-stage
CNLS-with-contextual-variables formulation of Johnson and Kuosmanen (2011).
The reported coefficient table, {cmd:e(b)}, and {cmd:e(V)} refer only to
these contextual-variable coefficients.  They do not provide sampling
uncertainty for the nonparametric frontier or for subsequent StoNED
efficiency scores.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{cmdab:in:puts(}{it:varlist}{cmd:)} specifies the input (regressor) variables.
At least one input variable is required.

{phang}
{cmdab:fr:ontier(}{it:string}{cmd:)} specifies the frontier orientation.
{opt frontier(production)} (the default) imposes concavity: the frontier
lies above the data.  {opt frontier(cost)} imposes convexity: the
frontier lies below the data.

{phang}
{opt context(varlist)} specifies contextual variables whose common linear
coefficients are estimated jointly with the CNLS frontier.  The option is
currently supported with {opt rts(vrs)}.  When {opt context()} is specified,
{cmd:convnls} reports an asymptotic coefficient table and stores {cmd:e(b)}
and {cmd:e(V)} for the contextual-variable coefficients.

{phang}
{opt vce(string)} specifies the variance estimator for coefficients reported
by {opt context()}.  {opt vce(ols)} is the default.  {opt vce(robust)} uses a
heteroskedasticity-robust sandwich variance matrix.

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
30 times this tolerance after two consecutive near-feasible cutting-plane
checks and the final feasibility solve.

{phang}
{opt rho(#)} sets the initial ADMM penalty parameter.  The default -1 lets
{cmd:convnls} choose the starting value automatically.  This is an advanced
diagnostic option and is not usually needed in routine use.

{phang}
{cmdab:maxi:ter(}{it:#}{cmd:)} sets the maximum number of ADMM inner
iterations per cutting-plane step.  Default is 20000.

{phang}
{cmdab:maxc:ut(}{it:#}{cmd:)} sets the maximum number of cutting-plane
outer iterations.  Default is 100.  Each cutting-plane iteration adds
violated Afriat constraints and re-solves the ADMM problem.  The algorithm
terminates early when no new violations are found.

{phang}
{opt prune(#)} removes constraints whose slack exceeds {it:#} before
each cutting-plane solve.  A value such as {opt prune(1e-4)} can speed
up computations for large samples by reducing the active constraint set.
Default is 0 (no pruning).

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
Euclidean distance d(i,h) in input space is within the p-th percentile of
distances from i are added to the initial set, where p is controlled by
{opt sweetpct()}.  An initial ADMM solve is performed on these constraints
before the cutting-plane loop starts, warm-starting the optimisation.
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
{opt vcstrat(all)} (the default) adds all currently violated constraints at
each iteration, which minimises the number of cutting-plane outer iterations
but may add a large number of constraints per step.
Use {opt addper(#)} to cap the number of newly added constraints to the
largest {it:#} violations per iteration.
{p_end}
{pmore}
{opt vcstrat(group)} implements the CNLS-G strategy of Lee et al. (2013):
for each observation i, only the most violated constraint (i, argmax_h violation)
is added per iteration.  At most n constraints are added per step.
Each step is faster but more iterations are needed.
This is recommended in combination with {opt init(sweet)} for large samples.
{p_end}
{pmore}
{opt vcstrat(one)} implements the CNLS-O strategy: only the single globally
most violated (i,h) pair is added per iteration.  This uses the fewest
constraints per step but requires the most iterations overall, and is
typically dominated by {opt vcstrat(group)}.
{p_end}

{phang}
{ul:Correctness guarantee.}
All combinations of {opt init()} and {opt vcstrat()} target the same
globally feasible CNLS objective solution as the default
{opt init(none) vcstrat(all)}, up to the stated numerical tolerances.
The cutting-plane loop terminates when no violated constraints remain,
irrespective of the initialisation or addition strategy (Lee et al. 2013,
Proposition 1).  Only runtime differs.
{p_end}


{marker remarks}{...}
{title:Remarks}

{pstd}
{ul:Algorithm.}
{cmd:convnls} uses a cutting-plane outer loop combined with an
Alternating Direction Method of Multipliers (ADMM) inner solver.
By default, the cutting-plane loop starts with no Afriat constraints and
iteratively adds all violated constraints until the solution is globally
feasible.
Constraint deduplication ensures each (i,h) pair is added at most once.
Optional CNLS+-G acceleration (Lee et al. 2013) is available through the
{opt init()} and {opt vcstrat()} options: {opt init(sweet)} pre-populates
the active set with nearby-observation pairs before the loop begins, and
{opt vcstrat(group)} adds only the worst violation per observation per step.
These options do not change the target CNLS problem; they reduce runtime for
n >= 100.

{pstd}
{ul:Scaling.}
The dependent variable and inputs are internally standardised
(mean zero, unit variance) before optimisation to improve numerical
conditioning.  Contextual variables specified in {opt context()} are centered
and standardised internally, and their reported coefficients are transformed
back to the original scale.  Fitted values and residuals are returned on the
original scale.

{pstd}
{ul:Saved variables.}
In addition to the variables named by {opt generate()} and {opt fitted()},
{cmd:convnls} saves {it:_convnls_alpha} (intercepts) and
{it:_convnls_beta1}, {it:_convnls_beta2}, ... (gradients) in the dataset.
These are overwritten if {cmd:convnls} is called again.

{pstd}
{ul:StoNED workflow.}
The standard two-stage StoNED estimator uses {help convnls} (Stage 1)
followed by {help stoned} (Stage 2):

{phang2}(1) Estimate CNLS and save residuals:{p_end}
{phang3}{cmd:convnls y, inputs(x1 x2) frontier(cost) generate(ehat)}{p_end}

{phang2}(2) Decompose residuals via method of moments and compute efficiency scores:{p_end}
{phang3}{cmd:stoned ehat, frontier(cost) generate(te)}{p_end}

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
Each ADMM inner solve is O(N^2 * m) per iteration.
For larger samples, estimation may take several minutes.
Use {opt verbose} to monitor progress.

{pstd}
{ul:Convergence warning.}
{cmd:e(converged)} is a status code: 0 means not converged, 1 means
converged to the specified tolerance, and 2 means practically converged with
a small remaining maximum violation.  The practical status uses the same
scale as {opt tolerance()}: the remaining maximum violation is below
30*{opt tolerance()} after repeated near-feasible checks and the final
feasibility solve.  If {cmd:e(converged) = 0}, consider
increasing {opt maxiter()} or loosening {opt tolerance()}.


{marker examples}{...}
{title:Examples}

{pstd}Copy this complete block to the Do-file Editor, omitting leading dots.
The examples generate their own data. Run from {cmd:preserve} through
{cmd:restore} to retain your original dataset on successful completion.
Estimation results and settings may change. Small samples demonstrate syntax;
they are not a recommended empirical research design.{p_end}
{* STONED_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}

{pstd}{ul:Example 1: Simulated cost frontier (VRS), StoNED two-stage}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. set obs 40}{p_end}
{phang2}{cmd:. gen ln_x1 = rnormal(1.5, 0.6)}{p_end}
{phang2}{cmd:. gen ln_x2 = rnormal(1.0, 0.4)}{p_end}
{phang2}{cmd:. gen u = abs(rnormal(0, 0.20))}{p_end}
{phang2}{cmd:. gen ln_c = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) + u}{p_end}

{phang2}{cmd:. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat) fitted(fhat)}{p_end}
{phang2}{cmd:. stoned ehat, frontier(cost) generate(te) inefficiency(uhat)}{p_end}
{phang2}{cmd:. summarize te uhat}{p_end}

{pstd}{ul:Example 2: Simulated production frontier (CRS)}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 4242}{p_end}
{phang2}{cmd:. set obs 100}{p_end}
{phang2}{cmd:. gen ln_x1 = rnormal(1.5, 0.6)}{p_end}
{phang2}{cmd:. gen ln_x2 = rnormal(1.0, 0.4)}{p_end}
{phang2}{cmd:. gen u = abs(rnormal(0, 0.20))}{p_end}
{phang2}{cmd:. gen ln_y = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) - u}{p_end}

{phang2}{cmd:. convnls ln_y, inputs(ln_x1 ln_x2) frontier(production) rts(crs) generate(ehat_p) init(sweet) sweetpct(3) vcstrat(group) maxcut(500)}{p_end}
{phang2}{cmd:. stoned ehat_p, frontier(production) generate(te_p)}{p_end}

{pstd}{ul:Example 3: Constant returns to scale}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. set obs 40}{p_end}
{phang2}{cmd:. gen ln_x1 = rnormal(1.5, 0.6)}{p_end}
{phang2}{cmd:. gen ln_x2 = rnormal(1.0, 0.4)}{p_end}
{phang2}{cmd:. gen u = abs(rnormal(0, 0.20))}{p_end}
{phang2}{cmd:. gen ln_c = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) + u}{p_end}
{phang2}{cmd:. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) rts(crs) generate(ehat_crs)}{p_end}

{pstd}{ul:Example 4: Verbose output and convergence diagnostics}

{phang2}{cmd:. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat_diag) verbose tolerance(1e-5)}{p_end}
{phang2}{cmd:. ereturn list}{p_end}

{pstd}{ul:Example 5: CNLS+-G accelerated algorithm (Lee et al. 2013)}

{phang2}{cmd:. * Recommended for n >= 100: Sweet Spot init + Group violation strategy}{p_end}
{phang2}{cmd:. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat_sweet) init(sweet) sweetpct(3) vcstrat(group)}{p_end}

{phang2}{cmd:. * Check which algorithm was used:}{p_end}
{phang2}{cmd:. ereturn list}{p_end}

{pstd}{ul:Example 6: Elementary Afriat init only}

{phang2}{cmd:. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat_ele) init(ele) vcstrat(group)}{p_end}

{pstd}{ul:Example 7: One-stage contextual-variable coefficient inference}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 270627}{p_end}
{phang2}{cmd:. set obs 60}{p_end}
{phang2}{cmd:. gen x = runiform()*10 + 1}{p_end}
{phang2}{cmd:. gen z = rnormal()}{p_end}
{phang2}{cmd:. gen u = abs(rnormal(0,.05))}{p_end}
{phang2}{cmd:. gen y = ln(x) + 2 + .35*z + rnormal(0,.05) - u}{p_end}

{phang2}{cmd:. convnls y, inputs(x) context(z) frontier(production) init(sweet) vcstrat(group) vce(robust) generate(ehat_ctx) fitted(fhat_ctx)}{p_end}
{phang2}{cmd:. ereturn list}{p_end}



{phang2}{cmd:. restore}{p_end}
{* STONED_HELP_EXAMPLES_END}{...}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:convnls} stores the following in {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(m)}}number of input variables{p_end}
{synopt:{cmd:e(q)}}number of contextual variables specified in {cmd:context()}{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}
{synopt:{cmd:e(sigma2)}}estimated variance (RSS / N){p_end}
{synopt:{cmd:e(converged)}}0 not converged, 1 converged, 2 practically converged{p_end}
{synopt:{cmd:e(cp_iters)}}number of cutting-plane iterations{p_end}
{synopt:{cmd:e(n_active)}}number of active Afriat constraints at convergence{p_end}
{synopt:{cmd:e(max_viol)}}maximum constraint violation at final solution{p_end}
{synopt:{cmd:e(addper)}}cap on newly added violations per iteration; 0 means no cap{p_end}
{synopt:{cmd:e(adamiter)}}Adam warm-start iterations{p_end}
{synopt:{cmd:e(adamlr)}}Adam warm-start learning rate{p_end}
{synopt:{cmd:e(adamlambda)}}Adam warm-start penalty weight{p_end}
{synopt:{cmd:e(adamscreen)}}Adam candidate-screening cap{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:convnls}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(inputs)}}names of input variables{p_end}
{synopt:{cmd:e(context)}}names of contextual variables, if specified{p_end}
{synopt:{cmd:e(frontier)}}{cmd:production} or {cmd:cost}{p_end}
{synopt:{cmd:e(rts)}}{cmd:vrs} or {cmd:crs}{p_end}
{synopt:{cmd:e(resid_var)}}name of the residual variable created{p_end}
{synopt:{cmd:e(fitted_var)}}name of the fitted-value variable created{p_end}
{synopt:{cmd:e(init_method)}}{cmd:none}, {cmd:ele}, {cmd:sweet}, or {cmd:adam}; initialisation used{p_end}
{synopt:{cmd:e(vcstrat)}}{cmd:all}, {cmd:group}, or {cmd:one}; violation strategy used{p_end}
{synopt:{cmd:e(vce)}}variance estimator used for {cmd:context()} coefficients{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}contextual-variable coefficient estimates, if {cmd:context()} is specified{p_end}
{synopt:{cmd:e(V)}}variance matrix for {cmd:e(b)}, if {cmd:context()} is specified{p_end}
{synopt:{cmd:e(V_robust)}}robust variance matrix for {cmd:e(b)}, if {cmd:context()} is specified{p_end}
{synopt:{cmd:e(cp_history)}}cutting-plane iteration history: iteration, active constraints, selected violations, newly added constraints or final-solve status, maximum violation{p_end}
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
bloom.rampike@gmail.com

{pstd}
Bug reports and suggestions are welcome.

{title:Software and article}

{pstd}Choonjoo Lee. {it:convnls and stoned: Stochastic Non-smooth Envelopment of Data in Stata}.
Journal information: {browse "https://www.stata-journal.com/":Stata Journal}.{p_end}
{pstd}Author-maintained software under the Academic Free License 3.0 (AFL-3.0).
See LICENSE and LICENSE_CODE.txt in the repository. Command versions and
subsequent help revisions are maintained separately.{p_end}
