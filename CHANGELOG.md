# Changes

## 2026-10-25 — Author update

- Updated the scalar commands and included their required internal helpers.
  The previous repository contained `convnls` 0.5.2 and `stoned` 0.3.2.
- **Changed `stoned` default:** `implementation(pystoned) method(mm) cet(addi)`
  reproduces the fixed pyStoNED 0.7.3 reference and requires `outcome()`.
  Existing scripts that intend the earlier Stata estimator family must add
  `implementation(legacy)`. Legacy scores use BC conditional expectations;
  public-reference scores use a different mapping. Numerical results can differ
  from older releases even with the explicit legacy mode.
- Added explicit parameter/score/optimizer diagnostics, residual provenance
  checks and output protection. `invalid(raw)` is a reproduction diagnostic,
  not a certification of an estimate. Read `help stoned` before using it.
- Context effects remain point estimates. Unvalidated covariance diagnostics
  are separated from ordinary inference; `e(V)` and confidence intervals are
  not supplied for these effects.
- Included the experimental DDF commands, optional sparse engine and CRS
  geometry diagnostics. DDF reference reproduction, directional legacy outputs
  and scalar 2012 methodology have distinct scopes. No new inference or
  Simar-Wilson estimator is introduced.

Computational versions remain in source headers. Help-only packaging edits
have a separate dated marker. Use the Git commit to identify the complete set.
