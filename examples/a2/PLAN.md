# A₂ / Boussinesq through sixth order

This document records the completed benchmark plan. The example pages
contain the implementation results and reproduction commands.

1. Fix the `κ=3` convention in [DZ (3.6.55)](https://arxiv.org/pdf/math/0108160#page=94)
   and the Boussinesq identification on [PDF p. 168](https://arxiv.org/pdf/math/0108160#page=168).
   Check the dispersive coefficient independently, either by the
   second cubic Gelfand–Dickey flow or from the DZ Poisson pencil.
2. In the flat variables `(V,U)` the exact differential equation is
   `U_t=V_x`, `V_t=3 U U_x+(ε²/4) U_xxx`. Its absence of higher input
   coefficients does not truncate the quasi-Miura map. Set
   `U=3b²`, `r,s=V±6b³`, so `r−s=12b³` and the diagonal speeds are
   `+3b,−3b`. Work locally where `b,r_x,s_x` do not vanish.
3. Use only that equation to compute the homogeneous Lie generators
   through ε⁶, via the coupled hodograph quadratures. The sparse
   Laurent state variable `b` keeps the cube-root speeds exact without
   replacing the physical canonical jets by differentiated powers of `b`.
   Verify both maps on the physical PDE through `O(ε⁷)`.
4. Recover the free energies from the dress map in the flat tau gauge:
   `U_(2g)=D_x² F_g`, `V_(2g)=D_x E_0 F_g`, and check both identities.
   Record the logarithmic primitive at genus one; no loop-equation
   coefficients are allowed as input to these integrations.
5. Compute the universal [DZ loop equation (3.10.95)](https://arxiv.org/pdf/math/0108160#page=150)
   independently from the A₂ Frobenius multiplication, intersecting
   metric, and periods. Substitute the computed densities and clear
   spectral denominators to test the coefficients at ε⁰, ε² and ε⁴.
   Verify symbolically for independent jets and a free spectral
   parameter, and include a negative control.
6. Archive the complete coefficients and a reproducible command. Keep all
   work finite-order and distinguish the exact input equation from its
   formal quasi-Miura output.
