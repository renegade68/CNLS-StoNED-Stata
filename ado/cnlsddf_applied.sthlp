{smcl}
{* Public packaging documentation revision 25oct2026; estimator unchanged.}{...}
{* *! applied-use guide ddf-interpretation1 25oct2026}{...}
{title:cnlsddf_applied -- Choosing directions and interpreting directional StoNED}

{p 4 4 2}
Applied-use companion to {help cnlsddf}, {help cnlsddf_sparse} and {help stonedddf}.
This guide adds no estimator options. The model is VRS with desirable outputs,
a common nonnegative direction, and scalar normal/half-normal directional errors.

{p 4 4 2}
The conditional u/v, signed distance D and targets discussed here require
{cmd:stonedddf, implementation(legacy)}. The default public-reference path returns
a coefficient-dependent ratio and a literal conditional quantity, using raw d
rather than the production residual -d. Its validity flags do not resolve this
sign convention or coefficient nonuniqueness. See {help stonedddf} before using
that reference output. It is not the directional target workflow in this guide.

{title:Choose a direction before inspecting efficiency results}

{p 4 4 2}
gx describes input contraction and gy output expansion. Use input-only directions
for contraction questions, output-only for expansion questions, or a joint
direction when both are relevant. A zero component fixes that coordinate during
the directional movement; it does not remove that input/output from estimation.
All components are common to observations and cannot be jointly zero.

{p 4 4 2}
Column means or medians can provide scales for a prespecified direction, but do
not establish its economic appropriateness. Calculate them on the intended
complete-case sample. When comparing subsamples, state whether the full-sample
direction is held fixed. Check that every proposed direction is nonzero.
Do not choose the direction that produces the most favorable efficiency or raw
log likelihood. Nonproportional directions change the distance and error model.

{p 4 4 2}
Match gx/gy order and units to inputs()/outputs(). A unit change in a quantity
requires the same change in its direction component. For g replaced by c*g,
c>0, distances, sigma_u, sigma_v and conditional u divide by c; original-coordinate
targets stay the same and log likelihood increases by N*log(c), up to numerical
tolerances. Directional u is neither a Farrell ratio nor a percentage.

{title:Read certification, boundary and domain results separately}

{p 4 4 2}
After stage one inspect e(converged). If it is not 1, stop before postestimation.
After stonedddf with implementation(legacy), save r(status), r(boundary), the scales and negative-coordinate
counts immediately: summarize, count and other r-class commands can overwrite r().
An error return certifies no new stage-two result; do not use stale r() or
silently switch methods. The example records rejected methods explicitly.

{p 4 4 2}
Interior (boundary=0) means the fitted decomposition has positive scales.
MM wrong_skew (1), zero_skew (2), legacy zero_variance (3; current stat1 rejects exact zero variance with r(498)) and likelihood normal_boundary
(4) are distinct outcomes. With sigma_u=0, all conditional u are zero: report the
boundary, and leave efficiency ranks or top groups undefined. These zeros do not
establish that all units are efficient. See {help stonedddf_legacy} for r(498) boundary
rejections and r(430) search exhaustion. Numerical success does not validate the
distributional assumptions, identify a causal effect or supply standard errors.

{p 4 4 2}
targets(t) gives (x-D*gx,y+D*gy). noiseadjusted(n) gives
(x+E(v|z)*gx,y-E(v|z)*gy). Moving the latter by E(u|z) along (-gx,+gy)
reaches the target; moving observed data by E(u|z) alone generally does not.
Distance D can be negative, so a fitted target need not be a contraction/expansion
improvement from the observed point. t_negative/n_negative flag any negative
coordinate. Report both counts and proportions. Do not clip, silently exclude,
or present flagged points as feasible recommendations. Envelope membership
does not by itself establish economic feasibility.

{title:Example with data already loaded}

{p 4 4 2}
For x1 x2 and desirable outputs y1 y2, after choosing gx(1 1) and gy(1 1)
in their quantity units:

{phang2}{cmd:cnlsddf, inputs(x1 x2) outputs(y1 y2) gx(1 1) gy(1 1) distance(d) residual(r) coefficients(c)}{p_end}
{phang2}{cmd:assert e(converged)==1}{p_end}
{phang2}{cmd:stonedddf, implementation(legacy) generate(s) targets(t) noiseadjusted(n)}{p_end}
{phang2}{cmd:return list}{p_end}

{p 4 4 2}
Check and save the returned boundary and negative-coordinate diagnostics before
other r-class commands. The chosen direction and admissibility of the MM
variance estimates must be assessed for the application.

{title:Scope}

{p 4 4 2}
Inference, observation-specific directions, nonnegative-target projection and
Simar-Wilson estimation are not supplied. See {help stonedddf_legacy} for the
formulas and {help stonedddf} for the distinct public-reference outputs.
