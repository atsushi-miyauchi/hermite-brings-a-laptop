# Code for the Hermite-Coefficient Certification Framework 

The code is organized to mirror the manuscript's current section structure.
The Hermite-coefficient certification is performed once, independently of any application, 
and the application-specific programs consume the resulting exact rational coefficient data.

## 1. Section 3: generate Hermite coefficient bounds once

The generic generator is

```bash
python hermite_coefficients/generate.py \
  --ks 2-32 --M 24 --R 8 --prec 256 --digits 40

python hermite_coefficients/generate.py \
  --ks 2-16 --M 32 --R 8 --prec 256 --digits 40
```

It implements the coefficient-certification pipeline developed in Section 3
of the manuscript: reduction to one-dimensional integrals, normalization and
symmetry grouping, finite-window/tail enclosures, and exact assembly of the
coefficient bounds.  It writes reusable JSON files containing certified
rational lower and upper bounds on `c_{k,m}`.

The repository already contains the reusable coefficient files used by the
current certificates:

- `hermite_coefficients/res/coeffs_k2-32_M24.json`
- `hermite_coefficients/res/coeffs_k2-16_M32.json`

Generation requires `python-flint`.  Once these JSON files have been produced,
the application verifiers in Section 4 use exact rational arithmetic only and
do not require `python-flint`.

The generator treats `m=0` exactly as `c_{k,0}=1/k`; numerical integration is
performed only for `m>=1`, as described after the proof in Section 3.4.

## 2. Section 4: exact application-specific verifiers

These programs do **not** regenerate Hermite coefficients.  They load one of
the reusable JSON files above and perform only the exact, problem-specific
certificate check appearing in Section 4.

### Section 4.3: MaxAgree

```bash
python maxagree/verify.py
```

This uses `M=24`, the 4-spokes/5-spokes mixture weights `0.364/0.636`, and
checks the `0.7818` positive-edge inequality by exact Bernstein arithmetic.
The negative-edge condition at `rho=0` is checked exactly.

### Section 4.4: MaxAgree[K]

Example for `K=3`:

```bash
python maxagree-k/verify.py --K 3 --target-ratio 0.8151
```

The verifier contains the Table 2 mixtures and uses `M=32`.  It checks the
positive-edge polynomial inequality by exact Bernstein arithmetic and the
negative-edge endpoint condition by exact rational evaluation.

### Section 4.5: Max K-Cut

Example for `K=4`:

```bash
python max-k-cut/verify.py \
  --K 4 --lower 0.857487 --upper 0.857488
```

The lower and upper bounds on the tight approximation ratio are checked
together from the same coefficient table, using `M=32`.

### Section 4.6: Modularity maximization

```bash
python modularity/verify.py --target-error 0.379
```

The current `0.3790` additive-error certificate requires `M=32`.

The near-optimality witness discussed later in Section 4.6 is checked by

```bash
python modularity_limitation/verify.py --target-value 0.3785
```

This uses `M=24` and checks the common witness `x=4/5` for spokes
`k=1,...,32`.

To run all exact application and noise-stability checks from Sections 4 and 5:

```bash
python verify_all.py
```

## 3. Section 5.4.3: quadratic bound on noise stability

```bash
python noise_stability/certify_gaussian_maxagree_bound.py
```

This verifies the numerical inequalities in the proof of
`Lambda_{3/4}(p) <= (3/4)*p^2 + (453/1000)*p` for all `p` in `[0,1]`,
where `Lambda_rho(p)` is the probability that a pair of rho-correlated
standard Gaussians both lie below `Phi^{-1}(p)`. The bound is used in the
MaxAgree integrality-gap construction with limiting ratio at most `0.802`.

The script uses exact rational intervals, alternating Taylor series, and
Machin's formula for pi to certify the two tangent inequalities at
`r=56/65`, using the proof's Jensen lower bound. The analytic argument in
Section 5.4.3 extends these checks to the whole interval. It prints the
certified enclosures and `PASS`, and requires only the Python standard
library, with no Hermite coefficient files.

## 4. Section 6: edge-wise upper-bound certificates

The LP certificate generators also load the same reusable Hermite coefficient
JSON files rather than regenerating them.

For MaxAgree:

```bash
cd edgewise_ub_code/maxagree_ub
python generate_certificate.py \
  --kfj '1-32' --kgw '1-5' --N 100000 --M 24 \
  --coeff-json ../../hermite_coefficients/res/coeffs_k2-32_M24.json

python check_certificate.py \
  res/maxagree_ub_N100000_M24_fj1-32_gw1-5.json \
  --require-below 0.78187
```

For MaxAgree[K]:

```bash
cd edgewise_ub_code/maxagreek_ub
python generate_certificate.py \
  --K-max 16 --N 100000 --M 32 \
  --coeff-json ../../hermite_coefficients/res/coeffs_k2-16_M32.json
```

The stored Section 6 certificates reference the shared coefficient files by
relative path and SHA-256 hash.

## Dependency split

- Section 3 Hermite coefficient generation: `python-flint`.
- Section 4 exact application verifiers: Python standard library only.
- Section 5 noise-stability verifier: Python standard library only.
- Section 6 certificate generation: `numpy`, `scipy`; MaxAgree with hyperplane
  roundings also uses `python-flint` for certified hyperplane probability
  bounds.
- Section 6 exact checking: exact rational arithmetic; the MaxAgree checker
  needs `python-flint` only if the optional GW re-audit is requested.
