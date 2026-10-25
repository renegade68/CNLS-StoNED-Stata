{smcl}
{* Public packaging documentation revision 25oct2026; estimator unchanged.}{...}
{title:cnlsddf_domain — CRS directional geometry diagnostics (research candidate)}

{p 4 4 2}
Version 0.1.0-domain2. Separate research command; requires Stata 16 or later,
macOS/Unix, and an external Python interpreter with NumPy. No optimization is
performed by this command. Use after certified sparse CRS crs2/dir1 estimation
with {cmd:coefficients()} specified.

{title:Syntax}

{p 8 8 2}
{cmd:cnlsddf_domain, python(}{it:executable}{cmd:)}
[{cmd:generate(}{it:prefix}{cmd:)} {cmd:replace}]

{p 4 4 2}
Default prefix: {cmd:_ddfdom}. New variables are double precision and missing
outside e(sample). All resulting names must fit Stata's 32-character limit.
{cmd:replace} allows exact output-name replacement, but never replacement of
bound first-stage variables. e() and the current data order are preserved.
Output failure restores the original data, types, labels and characteristics.

{title:Required source}

{p 4 4 2}
The current results must be from cnlsddf_sparse with rts(crs),
e(converged)==1, backend 0.1.0-crs2, bridge 0.1.0-dir1 and the CRS row-signature
schema. All input/output quantities, distances, residuals and coefficient
variables must remain available. The command verifies their joint-row binding
and direction. Whole-row sorting is allowed; changes outside the estimation
sample and unrelated variables do not invalidate the source.
The signature detects ordinary changes; it is not authentication of hostile
metadata or an independent optimizer certificate.

{title:Generated suffixes}

{p 4 4 2}
{cmd:_rawdist _envelope _zero}: saved distance, stored-plane envelope,
and ray zero using actual plane normalization. The original distance variable
is unchanged.

{p 4 4 2}
{cmd:_lower _upper}: bounds imposed by nonnegative quantities.
{cmd:_lbinf _ubinf}: lower bound is negative infinity / upper bound is positive
infinity. Unbounded limits are represented as missing plus these flags.

{p 4 4 2}
{cmd:_step _feasible _improve _observed}: maximal signed feasible step,
nonempty ray intersection, existence of a nonnegative step, and exact envelope
membership of the observed point. An empty ray has a missing step.

{p 4 4 2}
{cmd:_status}: 0 empty, 1 zero surface, 2 input cap, 3 near boundary.
{cmd:_endpoint}: 0 none, 1 zero surface, 2 input cap, 3 simultaneous zero and cap.
Near-boundary status does not override exact feasibility or invent an endpoint.

{p 4 4 2}
{cmd:_negany _negmat _negx _negy}: any raw negative quantity, materially negative
quantity, materially negative input, materially negative output.
Material threshold is -1e-10 after division by the saved-model quantity scales.
Scales are reconstructed as maximum absolute in-sample quantities, with zero
scales set to one, matching crs2.

{p 4 4 2}
{cmd:_minraw _marglo _marghi _minqty _minplane}: minimum scaled raw quantity,
scaled zero-to-lower and zero-to-upper margins, minimum scaled endpoint
quantity and endpoint plane margin. Infinite margins and absent endpoints
are missing; use bound and feasibility flags to distinguish them.

{p 4 4 2}
{cmd:_raw_x1} ... {cmd:_raw_y1} ... retain the original target quantities;
{cmd:_target_x1} ... {cmd:_target_y1} ... contain separate feasible endpoints.
Variable order follows e(inputs) and e(outputs). No quantity clipping or
replacement of signed distances is performed.

{title:Returned results}

{p 4 4 2}
r(N), r(empty), r(near_boundary), r(input_cap), r(raw_negative),
r(negative_step); r(cmd), r(version), r(source_status), r(generate).
r(empty) counts exact empty rays, including any near-boundary empty rays;
r(input_cap) counts status 2. Consequently these counts are not a mutually
exclusive partition. r(raw_negative) uses the material threshold.

{title:Interpretation}

{p 4 4 2}
A certified fit can have an empty physical ray. A feasible signed step may be
negative and is not a nonnegative improvement. Directional maximality does not
prove strong Pareto efficiency. Diagnostics concern the stored polyhedral CRS
technology, not a stochastic inefficiency decomposition. CRS MM/PSL/QLE
postestimation remains unsupported.

{title:Example after fitting}

{phang2}{cmd:cnlsddf_domain, python("/absolute/path/.venv/bin/python") generate(dom)}{p_end}
{phang2}{cmd:tabulate dom_status}{p_end}
{phang2}{cmd:list dom_rawdist dom_step dom_feasible dom_negmat if dom_feasible==0}{p_end}

{p 4 4 2}
First fit {help cnlsddf_sparse} with {cmd:rts(crs) coefficients(c)} and verify
{cmd:e(converged)==1}. Keep the fitted data and saved coefficients unchanged
before running this diagnostic. It does not provide CRS StoNED postestimation.
