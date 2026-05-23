# Alexander Invariants

A SageMath module for computing Alexander invariants of knots from PD or DT codes.
The module computes the Alexander module, its primary decomposition, invariant factors,
and the sequence of Alexander polynomials Δ₁(t), Δ₂(t), ... as introduced by Alexander.
Input may be given as a PD code, a DT code, or an alphabetical DT code.

---

## Definitions

The rational homology of the infinite cyclic cover of a knot is a module over
R = Q[t, t⁻¹]. It has an invariant factor decomposition as a direct sum of terms
R/<f_i(t)> where each successive f_i divides the previous one. The f_i can be
assumed to be primitive integer polynomials taking on a positive value at 0.

The n-th Alexander polynomial, or invariant, is the product of the f_i, starting
at f_n. The classical Alexander polynomial Δ(t) is thus Δ₁(t).

The primary decomposition splits the direct sum into primary pieces.

---

## Repository contents

| File | Description |
|------|-------------|
| `alexander_invariants.sage` | The main module. |
| `alexander_invariants_examples.ipynb` | A Jupyter notebook with worked examples. |
| `PD16Complete.txt.tar.bz2` | Names and PD notation for prime knots through 16 crossings. |
| `alexpolys_3_16_crossings.txt.tar.bz2` | Prime knots through 16 crossings, listed by name, rank of the rational Alexander module, and Alexander polynomials. |

---

## Dependencies

- [SageMath](https://www.sagemath.org/)
- [SnapPy](https://snappy.math.uic.edu/) — required only for DT input notation

---

## Loading

```sage
sage: load('alexander_invariants.sage')
```

SnapPy is imported automatically by the module. If you do not need DT input, SnapPy
is not required and the rest of the module will work without it.

---

## Main Function

```
alexander_invariants(code, decomp='alex_polys_decomp', notation='PD')
```

**Parameters**

| Parameter | Values | Default |
|-----------|--------|---------|
| `code` | PD code, DT code, or alphabetical DT code | — |
| `decomp` | `'alex_polys_decomp'`, `'primary_decomp'`, `'invariant_factor_decomp'` | `'alex_polys_decomp'` |
| `notation` | `'PD'`, `'DT'`, `'DT_alphabetical'` | `'PD'` |

---

## Output formats

All three decomposition formats encode the same underlying Alexander module;
they differ in how that information is presented.

### `'alex_polys_decomp'` (default)

Returns the list of Alexander polynomials [Δ₁(t), Δ₂(t), ...], each given
as a factored polynomial.  This is the format most directly comparable to
classical tabulations.

### `'primary_decomp'`

Returns the primary decomposition of the Alexander module. Each entry is of the form
`[f(t), [e₁, e₂, ...]]`, indicating cyclic summands Z[t,t⁻¹]/(fᵉⁱ).

### `'invariant_factor_decomp'`

Returns the invariant factors of the Alexander module as factored polynomials.

---

## Examples

### Basic example — 8₁₈, the first prime knot with noncyclic Alexander module

```sage
sage: pd = [[6,2,7,1],[8,3,9,4],[16,11,1,12],[2,14,3,13],
....:       [4,15,5,16],[10,6,11,5],[12,7,13,8],[14,10,15,9]]

sage: alexander_invariants(pd)
[(t^2 - 3*t + 1) * (t^2 - t + 1)^2, t^2 - t + 1]

sage: alexander_invariants(pd, decomp='primary_decomp')
[[t^2 - 3*t + 1, [1]], [t^2 - t + 1, [1, 1]]]

sage: alexander_invariants(pd, decomp='invariant_factor_decomp')
[(t^2 - 3*t + 1) * (t^2 - t + 1), t^2 - t + 1]
```
The primary decomposition `[[t^2 - 3*t + 1, [1]], [t^2 - t + 1, [1, 1]]]`
indicates that the Alexander module is

R/(t²−3t+1) ⊕ R/(t²−t+1) ⊕ R/(t²−t+1),

where R = Q[t, t⁻¹].

### DT notation input

```sage
sage: alexander_invariants([4, 8, -12, 2, -14, -16, -6, -10], notation='DT')
[(t^2 - t + 1) * (t^4 - t^2 + 1)]

sage: alexander_invariants('bdFaGHCE', notation='DT_alphabetical')
[(t^2 - t + 1) * (t^4 - t^2 + 1)]
```

Utility functions for converting between notations are also available directly:

```sage
sage: PD_from_dt([4, 8, -12, 2, -14, -16, -6, -10])
sage: PD_from_dt_alphabetical('bdFaGHCE')
```

### Connected sums

The module handles connected sums, for which the Alexander module is generally
non-cyclic. The following builds the 7-fold connected sum of the trefoil.

```sage
sage: pd = [[1,5,2,4],[3,1,4,6],[5,3,6,2]]    # trefoil
sage: k  = snappy.Link(pd)
sage: k2 = k.connected_sum(k)
sage: k4 = k2.connected_sum(k2)
sage: k6 = k4.connected_sum(k2)
sage: k7 = k6.connected_sum(k)

sage: print(factor(k7.alexander_polynomial()))   # classical single polynomial
(t^2 - t + 1)^7

sage: alexander_invariants(k7.PD_code())         # full module structure
[(t^2 - t + 1)^7, (t^2 - t + 1)^6, (t^2 - t + 1)^5, (t^2 - t + 1)^4,
 (t^2 - t + 1)^3, (t^2 - t + 1)^2, t^2 - t + 1]
```

### Unknot / trivial Alexander module

```sage
sage: pd = [[4,1,5,2],[3,1,4,6],[5,3,6,2]]
sage: alexander_invariants(pd, decomp='alex_polys_decomp')
[1]
sage: alexander_invariants(pd, decomp='primary_decomp')
[[1, [1]]]
sage: alexander_invariants(pd, decomp='invariant_factor_decomp')
[1]
```

---

## Return types

| `decomp` | Return type | Element type |
|----------|-------------|--------------|
| `'alex_polys_decomp'` | `list` | `Factorization` in `ZZ[t]` |
| `'primary_decomp'` | `list` of `[poly, list]` pairs | first entry: `Polynomial` in `ZZ[t]` |
| `'invariant_factor_decomp'` | `list` | `Factorization` in `ZZ[t]` |

---

## Input validation

The module checks that PD codes are well-formed and raises `ValueError` with a
descriptive message if not:

```sage
sage: bad = [[6,2,8,1],[8,3,9,4],[16,11,1,12],[2,14,3,13],
....:        [4,15,5,16],[10,6,11,5],[12,7,13,8],[14,10,15,9]]
sage: alexander_invariants(bad)
...
ValueError: Invalid PD code: crossing orientation condition failed at entry [6, 2, 8, 1]
```
