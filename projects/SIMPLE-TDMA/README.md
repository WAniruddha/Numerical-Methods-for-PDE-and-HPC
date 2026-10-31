# SIMPLE–TDMA: 2D Pressure-Correction and Scalar Transport (Fortran)

**Status:** Archived numerical-methods coursework/code; source files and outputs are retained separately in the original `SIMPLE_TDMA.zip` archive and are not included in this GitHub entry.

The archive contains Fortran 90 variants of a two-dimensional incompressible-flow solver using SIMPLE pressure–velocity coupling, a tridiagonal-matrix algorithm (TDMA), and passive-scalar/solute transport. The program reads domain dimensions, Reynolds and Péclet numbers, injection interval, Cartesian grid resolution, iteration limit, convergence tolerance, and under-relaxation coefficients from `parameters.par`.

**Evidence in archive:** `solver_SIMPLE-TDMA_ChannelSolute_v3.f90` and variants ending `_central`, `_central_press`, and `_central_pcorr`; example parameters (`Re=260`, `Pe=10`, `40 × 40` grid, `5000` maximum iterations, convergence criterion `1e-6`); a set of `.dat` output files and a TDMA lab PDF.

**Keywords:** Fortran 90, incompressible Navier–Stokes, finite-volume CFD, SIMPLE, TDMA, staggered-grid arrangements, pressure correction, scalar transport, residuals, under-relaxation.

**Disclosure:** Existing binaries/results are not supplied here; numerical accuracy, convergence, and validation have not been independently verified. Do not present specific validation outcomes without examining the outputs.