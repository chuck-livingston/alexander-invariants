# Alexander Invariants

This repository contains tools and data for working with Alexander invariants of knots: 
"alexander_invariants.ipynb" gives a basic tutorial. 
---

## Contact

For questions, comments, or collaboration, please contact:  
Charles Livingston – livingst@iu.edu

---

## Files

1. **Alexander_Invariants.pdf**  
   A paper that reviews the definition of Alexander invariants  
   and provides the theoretical background for the programs.

2. **alexander_invariants.sage**  
   The main program in Sage.

3. **alexander_invariants_examples.ipynb**  
   A Jupyter notebook with examples demonstrating the  
   use of the program, including optional parameters.

4. **pd_database.sqlite.zip**  
   A database for supplementary functions to find data  
   about prime knots of 16 crossings or less.  MUST BE  
   UNZIPPED FOR ACCESS TO 16 CROSSING DATABASE.  

5. **PD16Complete.txt**  
   A list of prime knots of 16 crossing or less with  
   names and indices.  MUST BE UNZIPPED.

---

## Basic Examples

```sage
sage: load('alexander_invariants.sage')

sage: pd = [[6,2,7,1],[8,3,9,4],[16,11,1,12],
...         [2,14,3,13],[4,15,5,16],[10,6,11,5],
...         [12,7,13,8],[14,10,15,9]]

sage: alexander_invariants(pd)
[(t^2 - 3*t + 1) * (t^2 - t + 1)^2, t^2 - t + 1]

sage: alexander_invariants(pd, decomp='primary_decomp')
[[t^2 - 3*t + 1, [1]], [t^2 - t + 1, [1, 1]]]

sage: alexander_invariants(pd, decomp='invariant_factor_decomp')
[(t^2 - 3*t + 1) * (t^2 - t + 1), t^2 - t + 1]

sage: alexander_invariants(pd, decomp='alex_polys_decomp')
[(t^2 - 3*t + 1) * (t^2 - t + 1)^2, t^2 - t + 1]

sage: index = find_index_from_name('K8a1'); print(index)
15

sage: pd = find_pd_from_name('K8a18'); print(pd)
[[6, 2, 7, 1], [8, 3, 9, 4], [16, 11, 1, 12],
 [2, 14, 3, 13], [4, 15, 5, 16],
 [10, 6, 11, 5], [12, 7, 13, 8], [14, 10, 15, 9]]

sage: alexander_invariants(pd)
[(t^2 - 3*t + 1) * (t^2 - t + 1)^2, t^2 - t + 1]


For more detailed documentation, see alexander_invariants_examples.ipynb and the source code.