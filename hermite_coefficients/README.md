# Reusable Hermite coefficient generator

`generate.py` implements the application-independent coefficient-certification
pipeline developed in Section 3 of the manuscript.  In particular, it follows
the reduction and grouping of Sections 3.1--3.2 and the rigorous enclosure and
assembly procedure of Sections 3.3--3.5.

Example:

```bash
python generate.py --ks 2-32 --M 24 --R 8 --prec 256 --digits 40
```

For each requested `k>=2` and `m=0,...,M`, the output contains exact rational
lower and upper bounds on `c_{k,m}`.  Degree zero is inserted exactly as
`c_{k,0}=1/k`, in accordance with the observation at the end of Section 3.4.
For `m>=1`, the program performs validated ACB integration on `[-R,R]`,
evaluates the analytic signed-tail bound in Arb, outward-rounds the result to
rational endpoints, applies the exact interval-square lower and upper bounds,
and assembles the coefficient in rational arithmetic.

The output is deliberately application-independent.  The exact verifiers in
Section 4 and the upper-bound certificates in Section 6 consume these reusable
coefficient files without rerunning the numerical integration.

Use `--details-csv FILE` to store the per-`(k,m,a,lambda)` interval data.
