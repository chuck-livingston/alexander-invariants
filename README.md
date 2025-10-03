# Alexander Invariants

This repository contains SageMath tools for computing Alexander invariants of knots. The Jupyter notebook 
"alexander_invariants.ipynb" gives a basic tutorial.  Be sure your kernel is SageMath. 


---

## Contact

For questions, comments, or collaboration, please contact:  
Chuck Livingston – livingst@iu.edu

---

## Files

1. **Alexander_Invariants.pdf**  
   A paper that reviews the definition of Alexander invariants  
   and provides the theoretical background for the programs.

2. **alexander_invariants.sage**  
   The main program in Sage.

3. **alexander_invariants_examples.ipynb**  
   A Jupyter notebook with examples demonstrating the  
   use of the program, including optional parameters.  Should
   have kernel set to SageMath.

--------

## Basic Examples in Sage

```
sage
sage: load('alexander_invariants.sage')

sage: pd = [[6,2,7,1],[8,3,9,4],[16,11,1,12],
...         [2,14,3,13],[4,15,5,16],[10,6,11,5],
...         [12,7,13,8],[14,10,15,9]]

sage: alexander_invariants(pd)  # = alexander_invariants(pd,decomp='alex_polys_decomp')
[(t^2 - 3*t + 1) * (t^2 - t + 1)^2, t^2 - t + 1]

sage: type(alexander_invariants(pd)[0])
<class 'sage.structure.factorization.Factorization'>

sage: alexander_invariants(pd, decomp='primary_decomp')
[[t^2 - 3*t + 1, [1]], [t^2 - t + 1, [1, 1]]]

sage: type(alexander_invariants(pd, decomp='primary_decomp')[0][0])
<class 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>

sage: alexander_invariants(pd, decomp='invariant_factor_decomp')
[(t^2 - 3*t + 1) * (t^2 - t + 1), t^2 - t + 1]

sage: type(alexander_invariants(pd, decomp='invariant_factor_decomp')[0])
<class 'sage.structure.factorization.Factorization'>


For more detailed documentation, see alexander_invariants_examples.ipynb and the source code.
