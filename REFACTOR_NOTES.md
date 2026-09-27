# Refactor notes

## Architectural change

Before the refactor, each application directory contained its own copy of the
ACB/Arb Hermite coefficient machinery and generated a problem-specific
coefficient JSON.  The Section 6 code contained another copy.

After the refactor, the code mirrors the current manuscript structure:

1. `hermite_coefficients/generate.py` is the only Hermite coefficient generator.
   It implements the application-independent coefficient certification developed
   in Section 3.
2. Its output is indexed only by the requested spokes parameters `k` and
   truncation degree `M` (plus numerical-generation parameters such as `R`,
   precision, and outward-rounding digits).
3. Each Section 4 application verifier loads a reusable coefficient JSON and
   performs only the exact problem-specific check.
4. The Section 6 LP generators also load the same coefficient JSON instead of
   regenerating Hermite coefficients.
5. Section 6 certificates record the coefficient file by relative path and
   SHA-256 hash, so the exact checker can verify that the intended reusable
   coefficient data were used.

This separation matches the logical flow of the paper: Section 3 produces
certified rational bounds on the Hermite coefficients, Section 4 consumes
those bounds to prove the algorithmic guarantees, and Section 6 reuses the
same bounds in the rounding-family upper-bound certificates.

Exploratory MaxAgree[K] mixture-estimation code and its diagnostic CSV files
have been removed from this package.  They are not part of any rigorous
certificate and are not needed to reproduce the manuscript results.

## Reusable coefficient data included

The two included files are exact repackagings of the certified rational data in
the supplied code; their coefficient arrays were checked entry-by-entry to be
identical to the originals.

- `hermite_coefficients/res/coeffs_k2-32_M24.json`
- `hermite_coefficients/res/coeffs_k2-16_M32.json`

## Verification performed

The refactored code was run against all manuscript values supplied in the code:

- Section 4.3, MaxAgree: 0.7818.
- Section 4.4, MaxAgree[K]: K=3,...,16, all displayed lower ratios.
- Section 4.5, Max K-Cut: K=3,...,16, both displayed lower and upper bounds.
- Section 4.6, modularity maximization: 0.3790.
- Section 4.6, modularity near-optimality: 0.3785 at x=4/5 for spokes
  k=1,...,32.
- Section 6, MaxAgree broad-family and Swamy-family stored certificates.
- Section 6, MaxAgree[K] stored certificates for K=3,...,16.

All of these checks pass after the refactor.

The shared coefficient generator itself was syntax-checked in this environment,
but was not numerically rerun here because `python-flint` is not installed in
the execution environment.  Its implementation is the generic ACB/Arb routine
from the supplied code, moved into one shared module.

## Important M observation

The supplied exact verifier shows that the 0.3790 modularity certificate in
Section 4.6 does not pass with M=24, but does pass with M=32.  By contrast, the
modularity near-optimality witness in the same subsection uses M=24.  The
refactored defaults therefore preserve those actually verified settings.
