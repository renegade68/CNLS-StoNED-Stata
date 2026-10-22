# CNLS and StoNED for Stata

Stata software accompanying Choonjoo Lee's manuscript,
*convnls and stoned: Stochastic Non-smooth Envelopment of Data in Stata*.
Journal information: [Stata Journal](https://www.stata-journal.com/).

`convnls` fits a shape-constrained Convex Nonparametric Least Squares (CNLS)
regression. `stoned` decomposes its residuals into noise and inefficiency and
computes conditional efficiency scores. Both production and cost orientations
are supported; returns to scale are specified in the first-stage CNLS fit.

The initial computational baseline is the final submission package dated
July 7, 2026: `convnls` 0.5.2 and `stoned` 0.3.2. Contact details and help examples
are maintained separately. Final journal citation and DOI will be added when
confirmed; later author-maintained software revisions are identified by Git commit.

## Install and begin

Requires **Stata 16 or newer**. No external Python, R or solver installation
is required. Download or clone this repository and keep the four files in
`ado/` together. In Stata:

```stata
adopath ++ "/your/path/CNLS-StoNED-Stata/ado"
which convnls
which stoned
help convnls
help stoned
```

Replace the path with your local folder. `which` shows the copy Stata will use
if an older installation is also present. This version uses a local-folder
installation and does not provide `net install` metadata. Access to a private
repository requires permission.

## Commands and examples

| Stage | Command | Purpose | Files |
| --- | --- | --- | --- |
| 1 | `convnls` | CNLS with monotonicity and convexity/concavity; VRS/CRS; optional contextual variables | [code](ado/convnls.ado) · [help](ado/convnls.sthlp) |
| 2 | `stoned` | Method-of-moments or pseudo-likelihood decomposition; efficiency and conditional inefficiency | [code](ado/stoned.ado) · [help](ado/stoned.sthlp) |

Each help file includes a complete block of synthetic-data examples. Copy it
into the Do-file Editor, omit leading dots, and run from `preserve` through
`restore`. The original dataset is restored on successful completion; estimation
results and settings may change. No separate dataset or do-file is required.

The examples cover cost/production frontiers, VRS/CRS, MOM/PSL, initialization
strategies, and contextual-variable inference. They demonstrate syntax and do
not replace the article's empirical replication or a study-specific validation.

## Interpreting results

Inspect `e(converged)` and `e(max_viol)` immediately after `convnls`, before
`stoned` replaces the estimation results. CNLS convergence codes are 0 (not
converged), 1 (converged), and 2 (practically converged); see the help for their
precise interpretation. Check `e(psl_converged)` when using `method(psl)`.
A successful command return alone does not establish numerical adequacy.

Use the same frontier orientation in both stages. `stoned` reports efficiency
scores between zero and one; its `inefficiency()` output is a different quantity.
The commands create or replace their requested output variables, so choose
output names that do not contain data you need to keep. See each help file for
assumptions, diagnostics and stored results.

## Versions and feedback

Computational versions are recorded in `.ado` headers. Help revisions have a
separate dated marker; help-only changes do not change the computational version.
Use a Git commit to identify the exact file combination.

For a problem report, use [GitHub Issues](https://github.com/renegade68/CNLS-StoNED-Stata/issues).
Include the Git commit, command, Stata version, `which` output, relevant
convergence diagnostics, and a small reproducible example showing expected
and actual results. Use synthetic data when research data are confidential.

## Attribution and license

The software is licensed under the [Academic Free License 3.0](LICENSE), following
the author's DEA software licensing policy. See [the software notice](LICENSE_CODE.txt).
Preserve author and method references. Cite the manuscript for this implementation
and the method papers listed in the help files for the underlying estimators.

This repository contains the two commands and their help. Manuscripts, signed
agreements, bank data, simulation scripts, Python comparison tools, logs and ZIP
archives remain outside this code repository. The software license does not
cover article manuscripts or publisher PDFs.

Contact: Choonjoo Lee — bloom.rampike@gmail.com
