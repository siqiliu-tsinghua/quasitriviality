# Toda / CP¹ loop-equation check

## Verification plan

1. Derive the CP¹ specialization from the period form of the Dubrovin–Zhang
   loop equation (3.10.95), using the periods and Gram matrix of Example 3.10.27.
   Record any differences from the printed (3.10.102).
2. Load the densities already reconstructed from the Toda quasi-Miura map in
   `toda-epsilon6.json`. Use no fitting or re-solving for their coefficients.
   Differentiate in flat jet coordinates by the exact chain rule, and hold the
   spectral parameter fixed in every total (x) derivative.
3. Compare the coefficients of `epsilon^0`, `epsilon^2`, and `epsilon^4` in
   the loop equation, with exact arithmetic and independent symbolic jets.
   The loop residual must be zero as a rational function of the spectral
   parameter, not just numerically at one or more sample points. Keep a
   reproducible command and a negative control.
