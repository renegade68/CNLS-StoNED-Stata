{smcl}
{* Public packaging documentation revision 25oct2026; estimator unchanged.}{...}
{* Interface documentation review 29sep2026; estimator ado unchanged.}{...}
{* *! version 0.3.3-dev5  28sep2026 (PSL normal boundary)}{...}
{* Documentation revision 28sep2026: score definitions and residual units; ado unchanged.}{...}
{viewerjumpto "Syntax" "stoned_legacy##syntax"}{...}
{viewerjumpto "Description" "stoned_legacy##description"}{...}
{viewerjumpto "Options" "stoned_legacy##options"}{...}
{viewerjumpto "Remarks" "stoned_legacy##remarks"}{...}
{viewerjumpto "Examples" "stoned_legacy##examples"}{...}
{viewerjumpto "Stored results" "stoned_legacy##results"}{...}
{viewerjumpto "References" "stoned_legacy##references"}{...}
{viewerjumpto "Author" "stoned_legacy##author"}{...}

{pstd}
Retained independent estimator. In this integrated interface every call must
specify {cmd:implementation(legacy)}. Public reproduction and raw-output policy
are documented in {help stoned}; this page describes the legacy formulas only.

{title:Title}

{phang}
{bf:stoned} {hline 2} StoNED Stage 2: Method of Moments or Pseudo-Likelihood decomposition

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:stoned} {it:resid_var} {ifin}
{cmd:,} {cmd:implementation(legacy)} [
{cmdab:front:ier(}{it:string}{cmd:)}
{cmdab:meth:od(}{it:string}{cmd:)}
{cmdab:gen:erate(}{it:newvar}{cmd:)}
{cmdab:ineff:iciency(}{it:newvar}{cmd:)}
{opt verbose}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmdab:front:ier(}{it:string}{cmd:)}}frontier type: {bf:production} (default) or {bf:cost}{p_end}
{synopt:{cmdab:meth:od(}{it:string}{cmd:)}}estimation method: {bf:mm} (default) or {bf:psl}{p_end}
{synopt:{cmdab:gen:erate(}{it:newvar}{cmd:)}}name for efficiency score variable (default: {it:_stoned_te}){p_end}
{synopt:{cmdab:ineff:iciency(}{it:newvar}{cmd:)}}name for conditional inefficiency estimate E(u|z){p_end}
{synopt:{opt verbose}}display moment and intermediate calculations{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:stoned} implements Stage 2 of the Stochastic Non-smooth Envelopment
of Data (StoNED) estimator of {help stoned_legacy##kuosmanen2012:Kuosmanen and Kortelainen (2012)}.
It decomposes the CNLS residuals produced by {help convnls} into a symmetric noise
component (v) and a one-sided inefficiency component (u) using the method of
moments, then computes firm-level efficiency scores using the conditional
expectation formula of {help stoned_legacy##battese1988:Battese and Coelli (1988)}.
Optionally, it also stores the observation-level conditional inefficiency
estimate E(u_i | z_i) of {help stoned_legacy##jondrow1982:Jondrow et al. (1982)}.

{pstd}
The model treats CNLS residuals as deviations from the conditional mean.
The command does not subtract their sample mean before PSL or scoring. After
estimating sigma_u, {cmd:stoned} reconstructs the composed error used in
the conditional efficiency formulas as epsilon_i + E(u) for cost frontiers
and epsilon_i - E(u) for production frontiers, where
E(u) = sigma_u*sqrt(2/pi).

{pstd}
The returned score is E[exp(-u)|z] for both orientations. Its interpretation
as an efficiency ratio requires u in natural-log units from a suitable log-error
model. With an outcome in levels, it is a unit-dependent exponential index;
{opt inefficiency()} instead reports E(u|z) in the outcome's units.
See {help stoned_legacy##units:Residual units and interpretation} below.

{pstd}
{it:resid_var} is the variable containing CNLS residuals, typically created by
the {opt generate()} option of {cmd:convnls}.  The typical two-stage StoNED
workflow is:

{phang2}{cmd:convnls} {it:depvar}{cmd:,} {cmdab:in:puts(}{it:varlist}{cmd:)}
  [{opt frontier(cost)}] {cmdab:gen:erate(}{it:ehat}{cmd:)}{p_end}
{phang2}{cmd:stoned} {it:ehat}{cmd:,} {cmd:implementation(legacy)} [{opt frontier(cost)}]
  {cmdab:gen:erate(}{it:te}{cmd:)}{p_end}


{marker options}{...}
{title:Options}

{phang}
{cmdab:front:ier(}{it:string}{cmd:)} specifies the frontier orientation,
matching the option used in {cmd:convnls}.
{opt frontier(production)} (the default) treats negative residual skewness as
evidence of a production-frontier composed error z = v - u.
{opt frontier(cost)} treats positive residual skewness as evidence of a
cost-frontier composed error z = v + u.

{phang}
{cmdab:meth:od(}{it:string}{cmd:)} specifies the Stage 2 estimation method.
{opt method(mm)} (the default) uses the method of moments (Kuosmanen and
Kortelainen, 2012), which solves closed-form moment equations for sigma_u and
sigma_v.
{opt method(psl)} uses the pseudo-likelihood (PSL) approach of Fan, Li, and
Weersink (1996), which maximises the normal--half-normal log-likelihood applied
to the reconstructed composed error z_i.  MOM estimates are used as starting values.
PSL requires numerical optimisation. Both methods impose the normal--half-normal
decomposition; finite-sample accuracy and boundary behaviour depend on the data.

{phang}
{cmdab:gen:erate(}{it:newvar}{cmd:)} names the new variable to hold the
BC exponential scores E[exp(-u)|z]. Scores lie on [0,1]. At a selected
zero-inefficiency boundary the score is 1 by convention, which is not evidence
that every unit is truly efficient. The economic ratio interpretation depends
on the residual units described below.
Default name: {it:_stoned_te}.

{phang}
{cmdab:ineff:iciency(}{it:newvar}{cmd:)} names the new variable to hold the
conditional inefficiency estimate E(u_i | z_i) of Jondrow et al. (1982).
The estimate is on the same scale as the CNLS residual.

{phang}
{opt verbose} prints the second and third central moments of the residuals
used in the MM calculation.


{marker remarks}{...}
{title:Remarks}

{pstd}
{ul:Method of Moments.}
Write s=+1 for cost and s=-1 for production. The centered model residual is
epsilon_i = v_i + s*(u_i-E(u_i)); the composed error is z_i=v_i+s*u_i.
Their population second and third central moments coincide:

{pmore}
m_2 = sigma_v^2 + sigma_u^2(1 - 2/pi)

{pmore}
m_3 = (4/pi - 1) * sqrt(2/pi) * sigma_u^3     [cost: positive; production: negative]

{pstd}
In a sample of N selected nonmissing observations, let
q_2=mean((epsilon-mean(epsilon))^2) and q_3=mean((epsilon-mean(epsilon))^3).
The implementation uses Stata's sample variance and skewness:

{pmore}
m_2 = [N/(N-1)]*q_2,   m_3 = [N/(N-1)]^(3/2)*q_3.

{pstd}
Thus {cmd:e(m3)} is skewness times sample_sd^3, not the unadjusted average
of cubed deviations. These are the program's finite-sample conventions.
With the expected sign s*m_3>=0, the formulas are:

{pmore}
sigma_u = (s*m_3 / c)^(1/3),   c = (4/pi - 1) * sqrt(2/pi) approx 0.2180

{pmore}
sigma_v = sqrt(m_2 - sigma_u^2*(1 - 2/pi))

{pstd}
If s*m_3<0, sigma_u is set to zero and sigma_v=sqrt(m_2); the absolute
value formula is not applied to wrong-sign skewness. If m_3 is exactly zero,
the same normal boundary follows. An arbitrarily small expected-sign m_3
is not rounded to zero. Nonpositive implied noise variance is handled below.

{pstd}
{ul:Efficiency scores.}
Given (sigma_u, sigma_v), firm-level efficiency scores are computed as the
conditional expectation E[exp(-u_i) | z_i], using the formula of
{help stoned_legacy##battese1988:Battese and Coelli (1988)}:

{pmore}
TE_i = exp(-mu*_i + sigma*^2/2) * Phi(mu*_i/sigma* - sigma*) / Phi(mu*_i/sigma*)

{pstd}
where sigma* = sigma_v * sigma_u / sigma, sigma = sqrt(sigma_v^2 + sigma_u^2),
z_i = epsilon_i + E(u) for cost and z_i = epsilon_i - E(u) for production,
and mu*_i = (+/-) sigma_u^2 * z_i / sigma^2 (+for cost, - for production).

{pstd}
If {opt inefficiency()} is specified, {cmd:stoned} also stores the JLMS
conditional inefficiency estimate

{pmore}
E(u_i | z_i) = mu*_i + sigma* phi(mu*_i/sigma*) / Phi(mu*_i/sigma*)

{pstd}
where phi and Phi are the standard normal density and distribution functions.

{marker units}{...}
{pstd}
{ul:Residual units and interpretation.}
{cmd:stoned} does not log, standardize, or recenter its input for scoring.
With natural-log outcomes and a suitable error specification,
exp(-u) is the production efficiency factor, or the frontier-to-actual cost
factor holding noise fixed. The stored BC score averages that factor
conditionally. It is generally different from exp(-E(u|z)); by Jensen's
inequality E[exp(-u)|z]>=exp(-E(u|z)).

{pstd}
With outcomes in currency, quantities, or another level unit, u and E(u|z)
have that same unit. The automatically returned exponential index then depends
on the chosen numerical unit and is not a unit-invariant output or cost ratio.
No additive efficiency-ratio formula is implemented by this command.
For a positive rescaling epsilon'=a*epsilon, the intended decomposition and
conditional mean satisfy sigma'_u=a*sigma_u, sigma'_v=a*sigma_v, U'=a*U,
whereas the new exponential score is E[exp(-a*u)|z] and generally changes.
Do not standardize residuals and interpret the resulting score as the original ratio.

{pstd}
Multiplying a level outcome by a constant before taking its natural logarithm
adds a constant to the logged outcome, rather than multiplying the log outcome.
Residual invariance then requires a model allowing the corresponding intercept
adjustment. If logs are supplied to {cmd:convnls}, its shape and
{opt rts(crs)} restrictions apply to those supplied coordinates; they do not
automatically impose shape or degree-one homogeneity in the original levels,
nor implement the separate multiplicative CNLS optimisation problem.

{pstd}
MM central moments are invariant to adding a constant to all selected residuals,
but the conditional scores and the PSL objective use the supplied residual
levels. Their invariance to a mean shift is not implied. CRS restrictions,
first-stage approximation error, or a Stage 2 subset can leave a nonzero
residual mean. Check the model and Stage 1 sample; silently centering the
input would change the estimator and, for marked residuals, invalidate provenance.

{pstd}
{ul:Pseudo-Likelihood (PSL).}
The PSL approach of {help stoned_legacy##fan1996:Fan, Li, and Weersink (1996)} maximises
the log-likelihood of the normal--half-normal composite error using the
mean-shifted composed-error reconstruction z_i(sigma_u):

{pmore}
z_i(sigma_u) = epsilon_i + sigma_u*sqrt(2/pi)     [cost]

{pmore}
z_i(sigma_u) = epsilon_i - sigma_u*sqrt(2/pi)     [production]

{pmore}
log L(sigma_u, sigma_v) = n*log(2) - n*log(sigma)
    + sum_i log phi(z_i(sigma_u)/sigma)
    + sum_i log Phi(lambda*z_i(sigma_u)/sigma)     [cost]

{pstd}
where sigma = sqrt(sigma_u^2 + sigma_v^2), lambda = sigma_u/sigma_v,
phi is the standard normal density, and Phi is the standard normal CDF.
For production, use Phi(-lambda*z_i(sigma_u)/sigma) inside the last logarithm;
the log Phi term is still added.
The shift is recomputed for each candidate sigma_u during likelihood
evaluation.  The optimisation is parameterised as
(log sigma_u, log sigma_v) to ensure positivity, and MOM estimates
serve as starting values. A third start (sigma_u,sigma_v)=(1.5,0.1) in RMSE units
checks a lower-noise region. Optimisation uses residual RMSE units and converts the
estimated standard deviations back to the original units. Residuals are not
recentered. The normal-boundary candidate is evaluated explicitly as described
below. PSL and MM need not give similar results in a small sample.

{pstd}
The MOM and lower-noise starts use BFGS; the neutral start uses Nelder--Mead
with log-parameter simplex deltas (0.1,0.1), which avoids a derivative-based
stopping rule on a flat near-zero-u likelihood. All starts have at most 500
iterations and ptol/vtol=1e{c -}9. BFGS additionally uses nrtol=1e{c -}10.

{pstd}
{ul:Cross-package definitions.}
Core {cmd:stoned} supports only {opt method(mm)} and {opt method(psl)}.
Its two-parameter PSL is not the one-parameter QLE implementation in the
fixed JSS pyStoNED 0.7.3 distribution. That QLE estimates lambda and sets
sigma^2=mean(epsilon^2)/[1-2*lambda^2/(pi*(1+lambda^2))]. The public MOM
uses population moments; on common admissible interior cases, the present
MM standard deviations are sqrt(N/(N-1)) times those from population moments.
Wrong-skew and unsupported-boundary policies also require separate checks.

{pstd}
The same public distribution returns (y-U)/y for ADDI production and (y+U)/y
for ADDI cost; MULT returns exp(-U) and exp(U), respectively. These are different
quantities from this command's E[exp(-u)|z], even with identical variance
parameters. Match moment conventions, score definitions, orientation and units
before comparing software. This describes the fixed 0.7.3 source, not every
later pyStoNED version or the directional-distance commands.

{pstd}
{ul:Wrong-skewness problem.}
If the sample skewness of the CNLS residuals has the wrong sign
(e.g., negative for a cost frontier), sigma_u cannot be identified via MOM.
{cmd:stoned} reports a warning and sets sigma_u = 0, assigning all variation to
noise.  This can occur with small samples or when true inefficiency is near zero.
For {opt method(psl)}, the optimiser may still find a solution with sigma_u > 0
in this case; users should interpret such results with caution.

{pstd}
{ul:Variance boundaries and numerical failures.}
At least three nonmissing residuals and finite moments with positive variance
are required. Samples smaller than three return r(2001); zero variance or
unrepresentable moments return r(498). For MM, sigma_u > 0 with implied
sigma_v^2 <= 0 is an unsupported zero-noise boundary and returns r(498).
It is not converted to TE=1 and E(u|z)=0. Those scores are valid for the
separate sigma_u=0, sigma_v>0 normal boundary, including the existing MM
wrong-skewness convention. Negative implied noise variance is not truncated.
PSL may still find a positive-noise solution when MM is invalid.
Conditional scores use log probabilities and continued fractions in the tails.
A nonfinite score calculation returns r(498), restoring the dataset.
Finite TE may round to 0 or 1 at machine precision.

{pstd}
{ul:Near-zero variance, skewness, and rankings.}
A numerically certified first-stage fit does not guarantee stable MM
decomposition near zero skewness or variance. With a constant outcome, an exact
zero residual can cause the required variance check to reject, while a solver's
tiny nonzero residual can lead to a formal normal-boundary score of 1. Neither
result alone demonstrates economically meaningful noise or full efficiency.
The command has no absolute near-zero-variance cutoff beyond its existing
validity checks, and does not silently clamp small positive skewness.

{pstd}
If sensitivity is material, compare certified fits from tighter first-stage
settings on the same sample and in the same units; record the residual differences,
MM boundary status, variance components and TE changes. {cmd:convnls}'s
{opt tolerance()} is a stopping criterion, not an observation-level residual-error bound.
Tighter settings can require more time and do not guarantee convergence or stability.
Almost equal scores should be reported as near ties or rank ranges with an explicit
numerical tolerance, alongside score differences. A changed rank or Spearman
correlation alone can exaggerate differences within such ties. A numerical tie
rule is not a statistical significance test, and {cmd:stoned} does not compute
rank uncertainty or efficiency confidence intervals.

{pstd}
{ul:PSL convergence and fallback.}
Only normally terminated, converged optimisations with a finite objective are
accepted as interior candidates. Iteration-limited and other unconverged iterates
are not PSL estimates. If all three starts fail, a valid MM decomposition is used; otherwise the command
returns r(498). {cmd:e(method)} records the requested method and
{cmd:e(effective_method)} records the method actually used. Displayed method
and variable labels use the effective method. {cmd:e(psl_converged)} is 0 or 1;
status 2 is no longer returned. Status 1 includes an explicit normal-boundary
selection after at least one interior start converged. A finite normal-boundary
likelihood alone cannot turn three failed starts into PSL success. These are numerical
diagnostics, not a proof of statistical identification or global optimality.

{pstd}
{ul:PSL normal-boundary selection.}
At sigma_u=0, the likelihood is maximised at sigma_v=sqrt(mean(epsilon_i^2)),
using the raw second moment and divisor N. For the best converged interior
candidate, {cmd:e(psl_ll_gain)} is (NLL_normal-NLL_interior)/N. If this gain is
at most 1e{c -}8, the command selects exactly sigma_u=0 and the normal-boundary
sigma_v; TE=1 and U=0. Otherwise it retains the interior candidate. This is a
dimensionless, per-observation numerical tie rule, not a significance test or
evidence that all units are truly efficient. Positive rescaling of residuals
does not change the intended selection rule. Exponential scores themselves are
not invariant to rescaling the units of u.

{pstd}
{cmd:e(psl_boundary_selected)} records the choice. {cmd:e(psl_selection)} is
{cmd:normal_boundary}, {cmd:interior}, or {cmd:mm_fallback}. A selected boundary
retains {cmd:e(effective_method)}={cmd:psl} and sets
{cmd:e(score_status)}={cmd:normal_boundary}. An interior candidate with
sigma_u/sigma_v < 1e{c -}4 receives a separate near-boundary note; it is not
silently set to zero. The zero-noise boundary sigma_v=0 remains unsupported.
The command also evaluates its limiting likelihood over the admissible half-normal
support. If this limit improves mean log-likelihood by more than 1e{c -}8 over the
selected candidate, the command returns r(498), saves no scores, and restores the
dataset. This rejects a known better unsupported limit; it does not certify a
global optimum or rule out every other interior mode.

{pstd}
{ul:Comparison with frontier MLE.}
Stata's built-in {helpb frontier} command maximises the full log-likelihood
and is appropriate when applied directly to observed outcomes.  When applied
to CNLS residuals, the MLE can converge to a boundary solution (sigma_u = 0)
if first-stage residuals are used without the mean-shift reconstruction needed
for a composed normal--half-normal error.
MM remains the closed-form default. A successful MM or PSL calculation does not
by itself establish correct distributional assumptions or stable boundary inference.


{pstd}
{ul:First-stage provenance.}
For residuals produced by {cmd:convnls} 0.5.3-dev8 or later with schema 1,
{cmd:stoned} checks the recorded convergence, inner-solver status, frontier and
residual-value signature before replacing outputs. Uncertified residuals,
incomplete/unsupported records, modified residual values or observation counts,
and a conflicting {opt frontier()} return r(498). Both strict (1) and practical
(2) CNLS certification are accepted and reported separately. The default
frontier remains production; cost residuals require {opt frontier(cost)}.

{pstd}
Records survive .dta save/use, rename, full clonevar, and other estimation commands.
Row reordering and display-format changes are allowed. {cmd:if}/{cmd:in} may select
a subset after checking the whole residual variable; physically dropping or
adding rows changes the signature and is rejected. Subsampling does not refit
CNLS or recenter residuals. Its statistical suitability is the user's responsibility.

{pstd}
An unmarked external or older residual is still accepted, with an explicit
unverified-source note. If the preceding estimation is CNLS and names that same
residual, its known nonconvergence or frontier mismatch is rejected. Otherwise
its origin cannot be inferred. This legacy check does not verify the residual
values. {cmd:generate} and CSV copies lose provenance and are treated as external
input. No source status for such input certifies the first stage.

{pstd}
The signature checks the residual multiset only. It does not authenticate user
edits to characteristics, match residuals to individual rows of other variables,
or check the current original X/Y. It protects ordinary workflows, not against
intentional removal or falsification of provenance.

{pstd}
{ul:Estimation sample.}
On normal completion, {cmd:e(sample)} identifies the residual observations used
in this Stage 2 fit: {cmd:if}/{cmd:in} restrictions with all missing residuals,
including extended missing values, excluded. Its count equals {cmd:e(N)}, and
TE and optional U are missing outside that sample. It describes Stage 2, which
may be a subset of the original CNLS sample counted in {cmd:e(stage1_N)}.
To keep the first-stage mask, generate a copy of {cmd:e(sample)} before calling
{cmd:stoned}. This posting creates no parametric coefficient vector or covariance
matrix and does not add generic predict or inference support.

{pstd}
Both commands support row reordering and in-memory {cmd:estimates store}/{cmd:restore}
for these sample masks. This is not a promise that masks accompany estimates
saved to a separate .ster file or arbitrary changes to the dataset. Errors are
not new estimates; check the return code before using previous e() results.

{pstd}
{ul:Output safety.}
The names in {opt generate()} and {opt inefficiency()} must be distinct and
must not match the residual input. Existing exact-name output variables are
replaced on normal completion; variable abbreviations never select a variable
for replacement. Invalid inputs or a runtime error leave the dataset unchanged.
An error does not provide a new successful estimate. Data protection uses
{cmd:preserve} and may require additional memory or temporary disk space.

{marker examples}{...}
{title:Examples}

{pstd}
The simulated examples below are self-contained. Save any data in memory
before running them. Empirical bank data and manuscript replication scripts
are distributed separately.

{pstd}
The ln_* variables below represent natural-log coordinates. CRS/VRS comparisons
illustrate the command's restrictions in those coordinates, not a test of
degree-one homogeneity in the original quantities. See
{help stoned_legacy##units:Residual units and interpretation} before applying the same
commands to outcomes measured in levels.

{pstd}{ul:Example 1: Simulated cost frontier}

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

{pstd}{ul:Example 2: Production frontier}

{phang2}{stata clear:. clear}{p_end}
{phang2}{stata set seed 4242:. set seed 4242}{p_end}
{phang2}{stata set obs 100:. set obs 100}{p_end}
{phang2}{stata gen ln_x1 = rnormal(1.5, 0.6):. gen ln_x1 = rnormal(1.5, 0.6)}{p_end}
{phang2}{stata gen ln_x2 = rnormal(1.0, 0.4):. gen ln_x2 = rnormal(1.0, 0.4)}{p_end}
{phang2}{stata gen u = abs(rnormal(0, 0.20)):. gen u = abs(rnormal(0, 0.20))}{p_end}
{phang2}{stata gen ln_y = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) - u:. gen ln_y = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) - u}{p_end}
{phang2}{stata convnls ln_y, inputs(ln_x1 ln_x2) frontier(production) rts(crs) generate(ehat_p) init(sweet) sweetpct(3) vcstrat(group) maxcut(500):. convnls ln_y, inputs(ln_x1 ln_x2) frontier(production) rts(crs) generate(ehat_p) init(sweet) sweetpct(3) vcstrat(group) maxcut(500)}{p_end}
{phang2}{stata stoned ehat_p, implementation(legacy) frontier(production) generate(te_p):. stoned ehat_p, implementation(legacy) frontier(production) generate(te_p)}{p_end}

{pstd}{ul:Example 3: CRS vs VRS comparison}

{phang2}{stata clear:. clear}{p_end}
{phang2}{stata set seed 42:. set seed 42}{p_end}
{phang2}{stata set obs 40:. set obs 40}{p_end}
{phang2}{stata gen ln_x1 = rnormal(1.5, 0.6):. gen ln_x1 = rnormal(1.5, 0.6)}{p_end}
{phang2}{stata gen ln_x2 = rnormal(1.0, 0.4):. gen ln_x2 = rnormal(1.0, 0.4)}{p_end}
{phang2}{stata gen u = abs(rnormal(0, 0.20)):. gen u = abs(rnormal(0, 0.20))}{p_end}
{phang2}{stata gen ln_c = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) + u:. gen ln_c = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) + u}{p_end}
{phang2}{cmd:. * CRS: alpha_i = 0 imposed}{p_end}
{phang2}{stata convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) rts(crs) generate(ehat_crs):. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) rts(crs) generate(ehat_crs)}{p_end}
{phang2}{stata stoned ehat_crs, implementation(legacy) frontier(cost) generate(te_crs):. stoned ehat_crs, implementation(legacy) frontier(cost) generate(te_crs)}{p_end}
{phang2}{cmd:. * VRS: free intercepts alpha_i}{p_end}
{phang2}{stata convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) rts(vrs) generate(ehat_vrs):. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) rts(vrs) generate(ehat_vrs)}{p_end}
{phang2}{stata stoned ehat_vrs, implementation(legacy) frontier(cost) generate(te_vrs):. stoned ehat_vrs, implementation(legacy) frontier(cost) generate(te_vrs)}{p_end}
{phang2}{cmd:. * CRS RSS >= VRS RSS (CRS is more constrained)}{p_end}
{phang2}{stata spearman te_crs te_vrs:. spearman te_crs te_vrs}{p_end}

{pstd}{ul:Example 4: MM vs PSL comparison}

{pstd}
Compare the method-of-moments (default) and pseudo-likelihood estimates
on the same CNLS residuals.

{phang2}{stata clear:. clear}{p_end}
{phang2}{stata set seed 42:. set seed 42}{p_end}
{phang2}{stata set obs 40:. set obs 40}{p_end}
{phang2}{stata gen ln_x1 = rnormal(1.5, 0.6):. gen ln_x1 = rnormal(1.5, 0.6)}{p_end}
{phang2}{stata gen ln_x2 = rnormal(1.0, 0.4):. gen ln_x2 = rnormal(1.0, 0.4)}{p_end}
{phang2}{stata gen u = abs(rnormal(0, 0.20)):. gen u = abs(rnormal(0, 0.20))}{p_end}
{phang2}{stata gen ln_c = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) + u:. gen ln_c = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) + u}{p_end}
{phang2}{stata convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat_mmpsl):. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat_mmpsl)}{p_end}
{phang2}{cmd:. * Method of moments (default)}{p_end}
{phang2}{stata stoned ehat_mmpsl, implementation(legacy) frontier(cost) method(mm) generate(te_mm) inefficiency(uhat_mm):. stoned ehat_mmpsl, implementation(legacy) frontier(cost) method(mm) generate(te_mm) inefficiency(uhat_mm)}{p_end}
{phang2}{cmd:. * Pseudo-likelihood}{p_end}
{phang2}{stata stoned ehat_mmpsl, implementation(legacy) frontier(cost) method(psl) generate(te_psl):. stoned ehat_mmpsl, implementation(legacy) frontier(cost) method(psl) generate(te_psl)}{p_end}
{phang2}{cmd:. * Compare sigma estimates and efficiency scores}{p_end}
{phang2}{stata di "MM  sigma_u=" e(sigma_u_mm)  " sigma_v=" e(sigma_v_mm):. di "MM  sigma_u=" e(sigma_u_mm)  " sigma_v=" e(sigma_v_mm)}{p_end}
{phang2}{stata di "PSL sigma_u=" e(sigma_u_psl) " sigma_v=" e(sigma_v_psl):. di "PSL sigma_u=" e(sigma_u_psl) " sigma_v=" e(sigma_v_psl)}{p_end}
{phang2}{stata spearman te_mm te_psl:. spearman te_mm te_psl}{p_end}

{pstd}{ul:Example 5: Bank cost efficiency (empirical)}

{pstd}
The full empirical example includes sample selection, log transformations,
CNLS estimation, StoNED decomposition, and saved results:

{pstd}The bank-data replication script is not included in this code repository.{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:stoned} stores the following in {cmd:e()}, including the function {cmd:e(sample)}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sigma_u)}}estimate of sigma_u (active method){p_end}
{synopt:{cmd:e(sigma_v)}}estimate of sigma_v (active method){p_end}
{synopt:{cmd:e(lambda)}}ratio sigma_u / sigma_v, if sigma_v > 0{p_end}
{synopt:{cmd:e(eu)}}estimated E(u) = sigma_u*sqrt(2/pi) used in the residual shift{p_end}
{synopt:{cmd:e(m2)}}sample variance N*q_2/(N-1){p_end}
{synopt:{cmd:e(m3)}}skewness times sample_sd^3; [N/(N-1)]^(3/2)*q_3{p_end}
{synopt:{cmd:e(sigma_u_mm)}}MOM estimate of sigma_u (always stored){p_end}
{synopt:{cmd:e(sigma_v_mm)}}MOM estimate of sigma_v (may be missing or zero when MM is invalid){p_end}
{synopt:{cmd:e(sigma_u_psl)}}PSL estimate or MM fallback of sigma_u (if {opt method(psl)}){p_end}
{synopt:{cmd:e(sigma_v_psl)}}PSL estimate or MM fallback of sigma_v (if {opt method(psl)}){p_end}
{synopt:{cmd:e(psl_converged)}}PSL status: 0 failed (valid MM fallback), 1 converged (if {opt method(psl)}){p_end}
{synopt:{cmd:e(psl_boundary_selected)}}1 if the explicit zero-u candidate was selected; otherwise 0 (PSL only){p_end}
{synopt:{cmd:e(psl_ll_gain)}}Best converged interior candidate's mean log-likelihood gain over the normal boundary; missing if all starts failed{p_end}
{synopt:{cmd:e(psl_boundary_tol)}}Per-observation numerical tie tolerance, 1e{c -}8 (PSL only){p_end}
{synopt:{cmd:e(stage1_converged)}}recorded first-stage status; missing for external input{p_end}
{synopt:{cmd:e(stage1_N)}}original first-stage N, distinct from the Stage 2 subset N; missing for external input{p_end}
{synopt:{cmd:e(mm_valid)}}1 if MM provides a supported positive-noise decomposition; 0 otherwise{p_end}
{synopt:{cmd:e(te_mean)}}mean efficiency score{p_end}
{synopt:{cmd:e(te_min)}}minimum efficiency score{p_end}
{synopt:{cmd:e(te_max)}}maximum efficiency score{p_end}
{synopt:{cmd:e(u_mean)}}mean conditional inefficiency, if {opt inefficiency()} specified{p_end}
{synopt:{cmd:e(u_min)}}minimum conditional inefficiency, if {opt inefficiency()} specified{p_end}
{synopt:{cmd:e(u_max)}}maximum conditional inefficiency, if {opt inefficiency()} specified{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(psl_selection)}}interior, normal_boundary, or mm_fallback (PSL only){p_end}
{synopt:{cmd:e(cmd)}}{cmd:stoned}{p_end}
{synopt:{cmd:e(method)}}requested {cmd:mm} or {cmd:psl}{p_end}
{synopt:{cmd:e(effective_method)}}actual {cmd:mm} or {cmd:psl}{p_end}
{synopt:{cmd:e(score_status)}}{cmd:interior} or {cmd:normal_boundary}{p_end}
{synopt:{cmd:e(source_status)}}{cmd:cnls_strict}, {cmd:cnls_practical}, {cmd:legacy_unverified}, or {cmd:external_unverified}{p_end}
{synopt:{cmd:e(stage1_frontier)}}recorded first-stage frontier, if known{p_end}
{synopt:{cmd:e(stage1_rts)}}recorded first-stage RTS, if known{p_end}
{synopt:{cmd:e(stage1_version)}}recorded CNLS version, if available{p_end}
{synopt:{cmd:e(resid_var)}}name of the CNLS residual variable{p_end}
{synopt:{cmd:e(frontier)}}{cmd:production} or {cmd:cost}{p_end}
{synopt:{cmd:e(te_var)}}name of the efficiency score variable created{p_end}
{synopt:{cmd:e(u_var)}}name of the conditional inefficiency variable created, if specified{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}
Dai, S., Y.-H. Fang, C.-Y. Lee, and T. Kuosmanen. 2024.
pyStoNED: A Python Package for Convex Regression and Frontier Estimation.
{it:Journal of Statistical Software} 111(6): 1-43.
DOI: {browse "https://doi.org/10.18637/jss.v111.i06":10.18637/jss.v111.i06}.{p_end}

{marker battese1988}{...}
{phang}
Battese, G. E., and T. J. Coelli. 1988.
Prediction of firm-level technical efficiencies with a generalized
frontier production function and panel data.
{it:Journal of Econometrics} 38(3): 387-399.
DOI: {browse "https://doi.org/10.1016/0304-4076(88)90053-X":10.1016/0304-4076(88)90053-X}.{p_end}

{marker jondrow1982}{...}
{phang}
Jondrow, J., C. A. K. Lovell, I. S. Materov, and P. Schmidt. 1982.
On the estimation of technical inefficiency in the stochastic frontier
production function model.
{it:Journal of Econometrics} 19(2-3): 233-238.
DOI: {browse "https://doi.org/10.1016/0304-4076(82)90004-5":10.1016/0304-4076(82)90004-5}.{p_end}

{phang}
Kuosmanen, T. 2008.
Representation theorem for convex nonparametric least squares.
{it:Econometrics Journal} 11(2): 308-325.
DOI: {browse "https://doi.org/10.1111/j.1368-423X.2008.00239.x":10.1111/j.1368-423X.2008.00239.x}.{p_end}

{marker fan1996}{...}
{phang}
Fan, Y., Q. Li, and A. Weersink. 1996.
Semiparametric estimation of stochastic production frontier models.
{it:Journal of Business and Economic Statistics} 14(4): 460-468.
DOI: {browse "https://doi.org/10.1080/07350015.1996.10524675":10.1080/07350015.1996.10524675}.{p_end}

{marker kuosmanen2012}{...}
{phang}
Kuosmanen, T., and M. Kortelainen. 2012.
Stochastic non-smooth envelopment of data: semi-parametric
frontier estimation subject to shape constraints.
{it:Journal of Productivity Analysis} 38(1): 11-28.
DOI: {browse "https://doi.org/10.1007/s11123-010-0201-3":10.1007/s11123-010-0201-3}.{p_end}


{marker author}{...}
{title:Author}

{pstd}
Choonjoo Lee{break}
Department of Defense AI/Robotics, Korea National Defense University{break}
Nonsan, Republic of Korea{break}
sarang90@korea.kr
