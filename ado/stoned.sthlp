{smcl}
{* Help revision: 2026-10-22.2; computational version unchanged.}{...}
{* *! version 0.3.2  26jun2026 (help examples and references verified)}{...}
{viewerjumpto "Syntax" "stoned##syntax"}{...}
{viewerjumpto "Description" "stoned##description"}{...}
{viewerjumpto "Options" "stoned##options"}{...}
{viewerjumpto "Remarks" "stoned##remarks"}{...}
{viewerjumpto "Examples" "stoned##examples"}{...}
{viewerjumpto "Stored results" "stoned##results"}{...}
{viewerjumpto "References" "stoned##references"}{...}
{viewerjumpto "Author" "stoned##author"}{...}

{title:Title}

{phang}
{bf:stoned} {hline 2} StoNED Stage 2: Method of Moments or Pseudo-Likelihood decomposition

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:stoned} {it:resid_var} {ifin}
[{cmd:,}
{cmdab:fr:ontier(}{it:string}{cmd:)}
{cmdab:met:hod(}{it:string}{cmd:)}
{cmdab:gen:erate(}{it:newvar}{cmd:)}
{cmdab:ineff:iciency(}{it:newvar}{cmd:)}
{opt verbose}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmdab:fr:ontier(}{it:string}{cmd:)}}frontier type: {bf:production} (default) or {bf:cost}{p_end}
{synopt:{cmdab:met:hod(}{it:string}{cmd:)}}estimation method: {bf:mm} (default) or {bf:psl}{p_end}
{synopt:{cmdab:gen:erate(}{it:newvar}{cmd:)}}name for efficiency score variable (default: {it:_stoned_te}){p_end}
{synopt:{cmdab:ineff:iciency(}{it:newvar}{cmd:)}}name for conditional inefficiency estimate E(u|z){p_end}
{synopt:{opt verbose}}display moment and intermediate calculations{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:stoned} implements Stage 2 of the Stochastic Non-smooth Envelopment
of Data (StoNED) estimator of {help stoned##kuosmanen2012:Kuosmanen and Kortelainen (2012)}.
It decomposes the CNLS residuals produced by {help convnls} into a symmetric noise
component (v) and a one-sided inefficiency component (u) using the method of
moments, then computes firm-level efficiency scores using the conditional
expectation formula of {help stoned##battese1988:Battese and Coelli (1988)}.
Optionally, it also stores the observation-level conditional inefficiency
estimate E(u_i | z_i) of {help stoned##jondrow1982:Jondrow et al. (1982)}.

{pstd}
CNLS residuals are centered around the estimated conditional mean.  After
estimating sigma_u, {cmd:stoned} reconstructs the composed error used in
the conditional efficiency formulas as epsilon_i + E(u) for cost frontiers
and epsilon_i - E(u) for production frontiers, where
E(u) = sigma_u*sqrt(2/pi).

{pstd}
{it:resid_var} is the variable containing CNLS residuals, typically created by
the {opt generate()} option of {cmd:convnls}.  The typical two-stage StoNED
workflow is:

{phang2}{cmd:convnls} {it:depvar}{cmd:,} {cmdab:in:puts(}{it:varlist}{cmd:)}
  [{opt frontier(cost)}] {cmdab:gen:erate(}{it:ehat}{cmd:)}{p_end}
{phang2}{cmd:stoned} {it:ehat}{cmd:,} [{opt frontier(cost)}]
  {cmdab:gen:erate(}{it:te}{cmd:)}{p_end}


{marker options}{...}
{title:Options}

{phang}
{cmdab:fr:ontier(}{it:string}{cmd:)} specifies the frontier orientation,
matching the option used in {cmd:convnls}.
{opt frontier(production)} (the default) treats negative residual skewness as
evidence of a production-frontier composed error z = v - u.
{opt frontier(cost)} treats positive residual skewness as evidence of a
cost-frontier composed error z = v + u.

{phang}
{cmdab:met:hod(}{it:string}{cmd:)} specifies the Stage 2 estimation method.
{opt method(mm)} (the default) uses the method of moments (Kuosmanen and
Kortelainen, 2012), which solves closed-form moment equations for sigma_u and
sigma_v.
{opt method(psl)} uses the pseudo-likelihood (PSL) approach of Fan, Li, and
Weersink (1996), which maximises the normal--half-normal log-likelihood applied
to the reconstructed composed error z_i.  MOM estimates are used as starting values.
PSL is generally more efficient than MM but requires numerical optimisation.

{phang}
{cmdab:gen:erate(}{it:newvar}{cmd:)} names the new variable to hold the
efficiency scores.  Scores are stored on the [0,1] scale, where 1 denotes full efficiency.
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
Under the half-normal model, the CNLS composite residual
epsilon_i = v_i +/- u_i has second and third central moments:

{pmore}
m_2 = sigma_v^2 + sigma_u^2(1 - 2/pi)

{pmore}
m_3 = (4/pi - 1) * sqrt(2/pi) * sigma_u^3     [cost: positive; production: negative]

{pstd}
Solving for sigma_u and sigma_v:

{pmore}
sigma_u = (|m_3| / c)^(1/3),   c = (4/pi - 1) * sqrt(2/pi) approx 0.2180

{pmore}
sigma_v = sqrt(m_2 - sigma_u^2*(1 - 2/pi))

{pstd}
{ul:Efficiency scores.}
Given (sigma_u, sigma_v), firm-level efficiency scores are computed as the
conditional expectation E[exp(-u_i) | z_i], using the formula of
{help stoned##battese1988:Battese and Coelli (1988)}:

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

{pstd}
{ul:Pseudo-Likelihood (PSL).}
The PSL approach of {help stoned##fan1996:Fan, Li, and Weersink (1996)} maximises
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
For a production frontier the sign of the last term is reversed.
The shift is recomputed for each candidate sigma_u during likelihood
evaluation.  The optimisation is parameterised as
(log sigma_u, log sigma_v) to ensure positivity, and MOM estimates
serve as starting values.
PSL is asymptotically more efficient than MM but requires numerical optimisation;
for small samples the two methods give similar results.

{pstd}
{ul:Wrong-skewness problem.}
If the sample skewness of the CNLS residuals has the wrong sign
(e.g., negative for a cost frontier), sigma_u cannot be identified via MOM.
{cmd:stoned} reports a warning and sets sigma_u = 0, assigning all variation to
noise.  This can occur with small samples or when true inefficiency is near zero.
For {opt method(psl)}, the optimiser may still find a solution with sigma_u > 0
in this case; users should interpret such results with caution.

{pstd}
{ul:PSL boundary solution.}
When {opt method(psl)} is used, the PSL log-likelihood may attain its
maximum near the boundary sigma_u = 0 in small or weakly skewed samples.
{cmd:stoned} evaluates the likelihood on the reconstructed composed error
z_i, but the first-stage CNLS fit can still leave limited one-sided
variation for Stage 2 in some applications.
{cmd:stoned} detects this condition (sigma_u < 1e{c -}4) and prints a note
recommending that users compare the result with {opt method(mm)}.  The MOM
estimator is more robust in this two-stage context and is the default for this
reason.

{pstd}
{ul:Comparison with frontier MLE.}
Stata's built-in {helpb frontier} command maximises the full log-likelihood
and is appropriate when applied directly to observed outcomes.  When applied
to CNLS residuals, the MLE can converge to a boundary solution (sigma_u = 0)
if first-stage residuals are used without the mean-shift reconstruction needed
for a composed normal--half-normal error.
The MOM estimator is more robust in this two-stage context and is the approach
recommended by {help stoned##kuosmanen2012:Kuosmanen and Kortelainen (2012)}.


{marker examples}{...}
{title:Examples}

{pstd}Copy this complete block to the Do-file Editor, omitting leading dots.
The examples generate their own data. Run from {cmd:preserve} through
{cmd:restore} to retain your original dataset on successful completion.
Estimation results and settings may change. Small samples demonstrate syntax;
they are not a recommended empirical research design.{p_end}
{* STONED_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}

{pstd}{ul:Example 1: Simulated cost frontier}

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

{pstd}{ul:Example 2: Production frontier}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 4242}{p_end}
{phang2}{cmd:. set obs 100}{p_end}
{phang2}{cmd:. gen ln_x1 = rnormal(1.5, 0.6)}{p_end}
{phang2}{cmd:. gen ln_x2 = rnormal(1.0, 0.4)}{p_end}
{phang2}{cmd:. gen u = abs(rnormal(0, 0.20))}{p_end}
{phang2}{cmd:. gen ln_y = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) - u}{p_end}
{phang2}{cmd:. convnls ln_y, inputs(ln_x1 ln_x2) frontier(production) rts(crs) generate(ehat_p) init(sweet) sweetpct(3) vcstrat(group) maxcut(500)}{p_end}
{phang2}{cmd:. stoned ehat_p, frontier(production) generate(te_p)}{p_end}

{pstd}{ul:Example 3: CRS vs VRS comparison}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. set obs 40}{p_end}
{phang2}{cmd:. gen ln_x1 = rnormal(1.5, 0.6)}{p_end}
{phang2}{cmd:. gen ln_x2 = rnormal(1.0, 0.4)}{p_end}
{phang2}{cmd:. gen u = abs(rnormal(0, 0.20))}{p_end}
{phang2}{cmd:. gen ln_c = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) + u}{p_end}
{phang2}{cmd:. * CRS: alpha_i = 0 imposed}{p_end}
{phang2}{cmd:. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) rts(crs) generate(ehat_crs)}{p_end}
{phang2}{cmd:. stoned ehat_crs, frontier(cost) generate(te_crs)}{p_end}
{phang2}{cmd:. * VRS: free intercepts alpha_i}{p_end}
{phang2}{cmd:. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) rts(vrs) generate(ehat_vrs)}{p_end}
{phang2}{cmd:. stoned ehat_vrs, frontier(cost) generate(te_vrs)}{p_end}
{phang2}{cmd:. * CRS RSS >= VRS RSS (CRS is more constrained)}{p_end}
{phang2}{cmd:. spearman te_crs te_vrs}{p_end}

{pstd}{ul:Example 4: MM vs PSL comparison}

{pstd}
Compare the method-of-moments (default) and pseudo-likelihood estimates
on the same CNLS residuals.

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. set obs 40}{p_end}
{phang2}{cmd:. gen ln_x1 = rnormal(1.5, 0.6)}{p_end}
{phang2}{cmd:. gen ln_x2 = rnormal(1.0, 0.4)}{p_end}
{phang2}{cmd:. gen u = abs(rnormal(0, 0.20))}{p_end}
{phang2}{cmd:. gen ln_c = 0.6*ln_x1 + 0.4*ln_x2 + rnormal(0,0.10) + u}{p_end}
{phang2}{cmd:. convnls ln_c, inputs(ln_x1 ln_x2) frontier(cost) generate(ehat_mmpsl)}{p_end}
{phang2}{cmd:. * Method of moments (default)}{p_end}
{phang2}{cmd:. stoned ehat_mmpsl, frontier(cost) method(mm) generate(te_mm) inefficiency(uhat_mm)}{p_end}
{phang2}{cmd:. * Pseudo-likelihood}{p_end}
{phang2}{cmd:. stoned ehat_mmpsl, frontier(cost) method(psl) generate(te_psl)}{p_end}
{phang2}{cmd:. * Compare sigma estimates and efficiency scores}{p_end}
{phang2}{cmd:. di "MM  sigma_u=" e(sigma_u_mm)  " sigma_v=" e(sigma_v_mm)}{p_end}
{phang2}{cmd:. di "PSL sigma_u=" e(sigma_u_psl) " sigma_v=" e(sigma_v_psl)}{p_end}
{phang2}{cmd:. spearman te_mm te_psl}{p_end}


{phang2}{cmd:. restore}{p_end}
{* STONED_HELP_EXAMPLES_END}{...}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:stoned} stores the following in {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sigma_u)}}estimate of sigma_u (active method){p_end}
{synopt:{cmd:e(sigma_v)}}estimate of sigma_v (active method){p_end}
{synopt:{cmd:e(lambda)}}ratio sigma_u / sigma_v, if sigma_v > 0{p_end}
{synopt:{cmd:e(eu)}}estimated E(u) = sigma_u*sqrt(2/pi) used in the residual shift{p_end}
{synopt:{cmd:e(m2)}}second central moment of residuals{p_end}
{synopt:{cmd:e(m3)}}third central moment of residuals{p_end}
{synopt:{cmd:e(sigma_u_mm)}}MOM estimate of sigma_u (always stored){p_end}
{synopt:{cmd:e(sigma_v_mm)}}MOM estimate of sigma_v (always stored){p_end}
{synopt:{cmd:e(sigma_u_psl)}}PSL estimate of sigma_u (if {opt method(psl)}){p_end}
{synopt:{cmd:e(sigma_v_psl)}}PSL estimate of sigma_v (if {opt method(psl)}){p_end}
{synopt:{cmd:e(psl_converged)}}PSL optimisation status: 0 failed, 1 fully converged, 2 partially converged (if {opt method(psl)}){p_end}
{synopt:{cmd:e(te_mean)}}mean efficiency score{p_end}
{synopt:{cmd:e(te_min)}}minimum efficiency score{p_end}
{synopt:{cmd:e(te_max)}}maximum efficiency score{p_end}
{synopt:{cmd:e(u_mean)}}mean conditional inefficiency, if {opt inefficiency()} specified{p_end}
{synopt:{cmd:e(u_min)}}minimum conditional inefficiency, if {opt inefficiency()} specified{p_end}
{synopt:{cmd:e(u_max)}}maximum conditional inefficiency, if {opt inefficiency()} specified{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:stoned}{p_end}
{synopt:{cmd:e(method)}}{cmd:mm} or {cmd:psl}{p_end}
{synopt:{cmd:e(resid_var)}}name of the CNLS residual variable{p_end}
{synopt:{cmd:e(frontier)}}{cmd:production} or {cmd:cost}{p_end}
{synopt:{cmd:e(te_var)}}name of the efficiency score variable created{p_end}
{synopt:{cmd:e(u_var)}}name of the conditional inefficiency variable created, if specified{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

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
bloom.rampike@gmail.com

{title:Software and article}

{pstd}Lee, C. Forthcoming. {it:convnls and stoned: Stochastic Non-smooth Envelopment of Data in Stata}.
Accepted for publication in the {browse "https://www.stata-journal.com/":Stata Journal}.{p_end}
{pstd}Author-maintained software under the Academic Free License 3.0 (AFL-3.0).
See LICENSE and LICENSE_CODE.txt in the repository. Command versions and
subsequent help revisions are maintained separately. The acceptance notice
refers to the article; it does not extend the article's validation claims
to later software changes.{p_end}
