# eLam

A Fortran Lambert solver using finite thin rectangular contour integrals, with elliptic, parabolic, hyperbolic and multiple-revolution coverage. This repository also preserves the supplied Exact v2.2 serial and vectorized kernels as independent, unchanged baselines.

The extended solver is a **new implementation of the rectangle-ratio/Newton strategy**, not a byte-identical version of Exact v2.2. It adds stable conic coordinates and branch safeguards. It does not call the confluent contour solver, ivLam, Izzo or Gooding. It is currently slower than ivLam in the mixed-domain comparison; the supplied original's speed claims must not be applied to this extension.

## Build and test

Requires CMake 3.20+ and GNU Fortran. OpenMP is used only to distribute independent original-kernel batch cases; use `-DELAM_OPENMP=OFF` for a serial build.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/elam_example
```

`elam` is a standalone library with no external data files or initialization. `elam_original` is a separate library for the preserved kernels and their adapters. The generated namespaced build copies change only the original module names, to permit both variants in one executable.

## API

```fortran
use elam_solver, only: elam_solve
! Consistent distance/time units; mu uses distance**3/time**2.
call elam_solve(r1, r2, tof, mu, direction, revolutions, branch, v1, v2, status)
! For exactly antipodal endpoints, supply the transfer plane explicitly:
call elam_solve(r1, r2, tof, mu, direction, revolutions, branch, v1, v2, status, plane_normal)
```

- `direction=+1` follows `r1 x r2`; `-1` follows the opposite orientation. It is not a global-z prograde flag.
- `revolutions>=0` counts complete revolutions. `branch=-1/+1` selects the lower/higher universal-anomaly root for multi-revolution transfers. For zero revolutions either value selects the unique root.
- `plane_normal` is perpendicular to both endpoint vectors. Its orientation, multiplied by `direction`, defines the plane for exact antipodal geometry. Coincident and same-ray geometry remain nonunique/unsupported.
- `elam_one` is the normalized (`mu=1`) entry used by the comparison, with the same argument order as `contour_one`; it also returns the geometry-conditioning flag `half`.
- Outputs are initialized to NaN on entry. Always check `status` before using them. Finite values alone are not an accuracy certificate.

| Status | Meaning |
|---:|---|
| 0 | Converged solution |
| 1 | Coalesced multi-revolution roots at minimum time |
| 5 | Time below the requested multi-revolution minimum |
| -1 | Invalid inputs |
| -2 | Degenerate geometry or missing/invalid antipodal plane |
| -3 | Nonunique/unsupported collinear geometry with supplied plane |
| -4 | Numerical convergence failure |

`half=-1/0/1` denotes degenerate/ordinary/near-collinear geometry. The dimensional API can solve exactly antipodal cases with a supplied plane; the comparison's C/MEX interface has no plane argument and reports those cases as unsupported.

## Method

The original Exact v2.2 uses a lookup-table guess in universal anomaly, an upper-half thin rectangle with eight Gauss–Legendre top nodes and one midpoint on each side, and a Newton correction. It is restricted to zero-revolution elliptic transfers.

The extension retains finite rectangular quadrature and Newton corrections. For each selected branch, it constructs a finite rectangle around the starter and uses the ratio of the first and zeroth contour moments of `1/(T-t)` to propose a root. A bracket rejects out-of-branch proposals and safeguards later Newton corrections with bisection. There is no confluent/Householder Lambert root step or alternate-solver fallback. Analytic parabolic limits and already-converged starters need no quadrature.

Additional coordinates are elliptic half anomaly, distance from the upper elliptic pole, hyperbolic half anomaly, and `sqrt(y)` for very short short-way hyperbolas. The last avoids subtracting nearly equal values at the hyperbolic branch point. Long-way hyperbolic time is rearranged algebraically to avoid cancellation. A separate Halley minimum search selects feasible multi-revolution branches. Radial/transverse velocity reconstruction avoids the near-antipodal cancellation of Lagrange coefficients. Logarithmic unit scaling supports very large and very small dimensional inputs.

Geometry, time derivatives, starters and velocity formulas were adapted from Lamberts-Comparison's `contour_simple`; the root iteration was replaced with finite rectangle quadrature and safeguarded Newton. This provenance is deliberate and does not constitute a fallback to that solver.

## Preserved originals

`original/exact_ser.F90`, `exact_vect.F90`, `exact_v1.F90`, `README.txt` and `CHANGES.txt` were copied byte-for-byte from the supplied `Exactv2.2` directory on 2026-10-22. `original/SHA256.json` and the `original_hashes` CTest verify the import. Exact v1 is retained for provenance but is not an active comparison entry. The supplied README contains historical results, not results reproduced by this repository.

Call `exact_initialize()` once, serially, before using `exact_one` or `exact_batch` from `exact_adapter`. The boolean first argument selects the vectorized kernel. These wrappers validate inputs, map direction conventions and reject multi-revolution requests (`-6`). They do not extend the original kernels or repair inaccurate finite answers. Batched arbitrary 3D inputs are rotated into local transfer planes before the original shared-plane API is called, and velocities are rotated back. That work is included in comparison timing.

No license declaration was supplied with these Exact source files; this package does not invent one. No ivLam source, coefficient data, compiled binaries or external submodules are required or bundled.

See [validation](docs/VALIDATION.md) for measured results and limits. The repository is [Ossama-Abdelkhalik-Lab/eLam](https://github.com/Ossama-Abdelkhalik-Lab/eLam). It can be developed independently or checked out as the `solvers/elam` submodule of Lamberts-Comparison; both use the same Git history. Commit and push solver changes here first, then record the updated submodule commit in the comparison repository.
