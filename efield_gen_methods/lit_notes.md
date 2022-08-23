# overview
Here, we expand on the two candidate methods for the electric field (EF) generation step in flock-it:
1. **GNN** EF < multipole expansion < multipole prediction via [Equivariant Multipole GNN](https://github.com/rinikerlab/EquivariantMultipoleGNN)
2. **GC-DNN** EF < electrostatic potential map predicted via [graph-convolutional DNN](https://github.com/AstexUK/ESP_DNN)

Our decision is based on the ability of each to reproduce the electric field such that we obtian reasonable STREUSEL surfaces of the molecules. This will be validated using the benchmark cases presented in the original STREUSEL publication, Table 1.

| system |
| ------ |



# to-do
[ ] - Gaussian SP to generate ESP on cx1 for each of the validation cases (Table 1)
[ ] - STREUSEL surface calculation for each validation case (Table 1)
[ ] - GNN surface generation [1]
[ ] - GC-DNN surface generation [2]
[ ] - error calculation (i.e. how far are the NN-generated EF surfaces from the "true" EF surface (STREUSEL). This will be visualized in a manner similar to that presented in [1].


# 1. EF from multipole expansion
## identified problem
describing electrostatic interactions is challenging and the current, common method (fixed partial-charge approximation) fails at short range because it doesn't account for conformational changes or anisotropic effects
implementing a multipole expansion enables an exact treatment of the electrostatic potential.
obtaining an ESP requires QM calculations, which are expensive
this work presents an equivariant GNN to solve this problem.
the published model eliminates the need for QM calculations by predicting atomic multipoles (up to teh quadrupole)
the equivariant architecture enforces the correct symmetry without needing local reference frames
the published GNN reproduces ESP for several benchmark systems.

## theory
### electrostatic potential and multipoles
The ESP map of a charge distribution ($\rho(r')$) at a point ($r$) is defined as

$ Ves(r) = \int\frac[\rho(r')][|r - r'|]dV' $

The Taylor expansion of the ESP map around the center of the charge distribution yields the multipole expansion

$ Ves(r) = \frac[M^{(0)}][|r|] + \frac[M_a^{(1)}r_a][|r|^3] + \frac[M_{\alpha\beta}^{(2)}(3r_\alpha r_\beta - r^2\delta_{\alpha\beta}][2|r|^5] + \ellip $

where $\alpha$ and $\beta$ run over the coordinates of the system and the multipoles are given as

| term | pole | colloquial symbol | 
| ---- | ---- | ----------------- |
| $M^{(0)}$ | $q$ | monopole | 
| $M^{(1)}$ | $\mu$ | dipole |
| $M^{(2})$ | $\theta$ | quadrupole |

It should be noted that the prefactor $(4\pi\epsilon_0)^{-1}$ is omitted for clarity.

Since molecular charge densities can be decomposed into atomic charge densities,

$ \rho(r) = \sum_i \rho_i(r) $

we can apply the same concent to the atomic densities, which yields atomic (or distributed) multipoles.

# atomic multipoles
From several available methods to calculate atomic multipoles the authors elected to use the minimal basis iterative Stockholder (MBIS) method. This decision was made because MBIS features a low computational cost and fast convergence wrt to the multipole order.
The MBIS method falls under the atom-in-molecule method class, which partitions total charge densities into atomic contributions. The MBIS method is based on the minimization of the Kullback-Leibler (KL) divergence between a pro-density ($\rho^0_A(r)$) based on a minimal expansion in atom-centered s-type Slater functions and atarget molecular density ($\rho_A(r)$),

$ \rho^0_A(r) = \sum^{m_A}_{i=1}\rho_{Ai}^0(r) $

where the pro-density ($\rho_A^0(r)$) is the pro-atomic density constructed from the Slater functions,

$ \rho^0_{Ai}(r) = N_{Ai}f_{Ai}(r) = \frac[N_{Ai}][\sigma_{Ai}^38\pi]\exp(-\frac[|r-R_A|][\sigma_{Ai}]) $

where $N_{Ai}$ is a fitting parameter for the electron population ($\sigma_{Ai}.



# 2. EF from ESP map


# References
| Ref. No. | DOI |
| -------- | --- |
| 1 | 10.1021/acs.jctc.1c01021 |
| | |
| 2 | 10.1021/acs.jmedchem.9b01129 |
