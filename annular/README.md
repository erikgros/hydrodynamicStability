# I did this as part of the instabilities course project between: 26.11.16 - 14.01.17

* Y. Renardy (1997): "In the incompressibility condition, v and w need to be expanded to the same degree polynomial as u' so the expansions for v and w are truncated at the (Ni 1)th Chebyshev polynomial.In the radial component of the equation of motion, p' is required to the same accuracy as u' so p is also truncated at the (Ni 1)th polynomial. Together with , the total is 4N1 + 4N2 + 3 unknowns. In the core, we express the continuity equation, the radial component of the momentum equation and the other two components of the momentum equation, respectively, to degree N1, (N1 1) and (N1 2), while in the annulus, the momentum equations are expanded to one degree lower"
* real part of the energy equation does not have any viscous dissipation (D) nor any surface tension (B1), it only contains some I and some B2
* The parameters of Fig2 (paper II) seem to produce physically spurious eigenvalues (see Boyd, Chap7), there are eigenvalues of infinite magnitude irrespective of the resolution. The computed eigenvalues approach them by increasing with resolution, while the corresponding eigenfunctions are sawtooth oscillations.
* The boundary conditions at r=0 can not be ignored, since the governing equation has singular coefficients at r=0 it cannot be imposed there.
* Tab 1 values are approached by both sysI and sysII but sysI converges more robustly and to closer values.
* RUNNING WITH DIFFERENT RESOLUTIONS SEEMS TO BE A GOOD SOLUTION BUT EFFICIENT IMPLEMENTATION is not obvious!

## Lessons Learned:
* High order derivatives cause ill conditionning of Chebyshev pseudospectal matrices
* The vander(x) matrix is well conditioned over the Chebyshev nodes provided there are not too many. This limits the possibility of obtaining better results by increasing the number of collocation points.
* The only place that the density ratio Czenters into the equations is through the jump in the perturbation pressure in the normal stress balance equation at the interface??? (from Lub pipelining 4 chap 6)

## sysI
* It took ~4 days to implement with qualitative agreement and more then 2 weeks to debug!
* There are large spurious eigenvalues!
* it is noticeably slower than sysII
* cond(A) smaller as for sysII but cond(B)=Inf!!!!!
* There seems to bee a serious dependence on the number of points in the outer region. Reducing the number of points makes the eigenvalues seamingly switch to a different capillary like branches.
* Fixing the problem in the normal stress BC made almost no difference for the results.
* High Re (Rei~1e4for RP and Chandra works just as well. But it has an effect on Chandra where the curves are never perfectly matching the analytical dispersion relation.
* The dWzetaJump seems not to affect (even the zeta~1000 Why?) results!
* There is spurious peak at interface (between two layers of the same fluid) when trying to simulate single-phase pipe flow

## sysII (axi and for equal densities), NOT WORKING!
* took me almost 2 weeks to implement and get qualitative agreement
* There are large spurious eigenvalues !
* Condition numbers are very large!
