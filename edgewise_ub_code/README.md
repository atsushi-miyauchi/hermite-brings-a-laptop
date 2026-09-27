# Edge-wise upper-bound certificate code

These programs implement the finite-dual certification framework of Section 6.
Hermite coefficient generation is no longer duplicated here.  Both generators
load coefficient bounds produced by `../../hermite_coefficients/generate.py`.

- `maxagree_ub/`: MaxAgree with an arbitrary finite family of spokes/FJ and
  hyperplane/GW roundings.
- `maxagreek_ub/`: MaxAgree[K] with spokes `k=1,...,K`.

Certificate generation uses a floating-point LP only to locate a candidate dual
support.  The support is rationalized and normalized exactly, and the stored
certificate is checked independently in rational arithmetic.  Each certificate
records the SHA-256 and relative path of the reusable coefficient JSON on which
it depends.
