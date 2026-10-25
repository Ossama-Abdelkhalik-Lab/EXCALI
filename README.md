# eLam

A Fortran Lambert solver using finite thin rectangular contour integrals, with elliptic, parabolic, hyperbolic and multiple-revolution coverage. This repository also preserves the supplied Exact v2.2 serial and vectorized kernels as independent, unchanged baselines.

The extended solver is a **new implementation of Exact v2.2's thin-rectangle ratio**, not a byte-identical version of it. It takes the rectangle's zero-height limit, which needs only real arithmetic, applies it in one or two fixed stages with no Newton polish, and adds stable conic coordinates and branch safeguards. It does not call the confluent contour solver, ivLam, Izzo or Gooding. On the mixed benchmark it takes 0.4–1.1 times ivLam's time per transfer on Apple silicon (0.91 overall) and 0.6–1.4 times on a Linux x86-64 runner, at equal or better accuracy on both. It is portable standard Fortran; CI covers Linux x86-64, Linux ARM64 and macOS (see [validation](docs/VALIDATION.md)); the supplied original's speed claims must not be applied to this extension.

## Build and test

Requires CMake 3.20+ and GNU Fortran. OpenMP is used only to distribute independent original-kernel batch cases; use `-DELAM_OPENMP=OFF` for a serial build.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/elam_example
```

An optional profile-guided build makes the solver about 4% faster. It uses the same flags otherwise, including `-fno-fast-math`:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DELAM_PGO=GENERATE
cmake --build build -j && ./build/elam_pgo_train
cmake -S . -B build -DELAM_PGO=USE && cmake --build build -j
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

A full write-up, from the Negrete–Abdelkhalik paper through Exact v1 and v2.2 to this code, with a flowchart, the code breakdown and the measured results, is in [docs/method/eLam_method.pdf](docs/method/eLam_method.pdf). Its LaTeX source is alongside (`latexmk -pdf eLam_method.tex`).

The original Exact v2.2 uses a lookup-table guess in universal anomaly, an upper-half thin rectangle with eight Gauss–Legendre top nodes and one midpoint on each side, and a Newton correction. It is restricted to zero-revolution elliptic transfers.

The extension keeps the thin rectangle but evaluates its height-to-zero limit. On a real interval `[c-L, c+L]`, with `f = 1/(T-t)` and Gauss–Legendre nodes `x_i` and weights `w_i`,

```
root = c + Q[((x-c) f)'] / Q[f'],    Q[g'] = g(c+L) - g(c-L) - L * sum_i w_i g'(x_i)
```

`Q` is the quadrature error of an exact derivative. The ratio is exact for a simple pole anywhere in the interval, so only the analytic remainder of `1/(T-t)` limits it, as `(L/2R)^(2n)` for the nearest other singularity at distance `R`. It needs `T` at the nodes and endpoints and `T'` at the nodes; `T'` comes from the closed-form ODE derivative. There is no complex arithmetic. Node trigonometry uses angle addition with compile-time constants for power-of-two half-widths, and the half-angle functions of the centre are carried between stages, and into the velocities, by angle addition. Each coordinate has one fused node loop.

Each solve runs a fixed pipeline: closed-form starter, an 8-node stage with `L ~ R/4`, and, only when the stage's own error estimate is above rounding, a 3-node stage centred on the first estimate. The error estimate is the disagreement of the two independent ratios `Q[((x-c)f)']/Q[f']` and `Q[((x-c)^2 f)']/Q[f']`. `R` is computed analytically: branch points of `sqrt(y)`, the mirror root of the even zero-revolution time function, the other multi-revolution branch, and the complex roots of `T = t` beside a pole of `T` (at `sqrt(3)` times the pole distance). Every accepted interval must show a sign change of `T - t` without containing a pole. An estimate outside its interval re-centres, and an inconsistent one widens the interval. The rare cases whose starter still misses (at most 0.8% by regime) continue as a bracketed contour: the sign of `T - t` at the centre and the rectangle's own sign test shrink the starter's bracket, and every estimate is still a zero-height rectangle ratio. Very short hyperbolas (`t < 0.1 t_p`) use the same bracketed contour in `sqrt(y)`. The root solve has no Newton step and no complex arithmetic; the only derivative-based iteration is the one- or two-step Halley search for the multi-revolution minimum time. Analytic parabolic limits need no quadrature.

Coordinates are `X = ±x²` near the parabola (a power series, which removes the mirror root), elliptic half anomaly, distance from the upper elliptic pole, hyperbolic half anomaly, and `sqrt(y)` for very short short-way hyperbolas. The last avoids subtracting nearly equal values at the hyperbolic branch point. Long-way hyperbolic time is rearranged algebraically to avoid cancellation. For multiple revolutions, every ellipse through both points has `a >= s/2`, so `T > m*2*pi*(s/2)^1.5` rejects most infeasible times before `T` is evaluated. Otherwise a Halley search for the minimum time (one or two steps, from the large-`m` minimum of `y^1.5/sin^3 v`) stops as soon as its error bound decides feasibility and places the minimum well inside the distance to the root. Within 1% of the minimum time, the branch starter inverts the local cubic Taylor model. Starters need only a few digits, so they use a bit-pattern cube root with one Halley correction. Branch-point distances use a fast logarithm. Near a pole of `T`, starters solve the pole asymptote `T ~ n*pi*y^1.5/(sqrt2*d^3)` with `y` quadratic in the pole distance `d`, which also covers the plateau of near-full-period long-way transfers. Radial/transverse velocity reconstruction avoids the near-antipodal cancellation of Lagrange coefficients. Logarithmic unit scaling supports very large and very small dimensional inputs.

Geometry, time derivatives, starters and velocity formulas were adapted from `contour_simple`, a short variant of ConfluentContour that has since been removed (it remains in that repository's history); the root iteration was replaced with zero-height rectangle quadrature. eLam does not call ConfluentContour.

## Preserved originals

`original/exact_ser.F90`, `exact_vect.F90`, `exact_v1.F90`, `README.txt` and `CHANGES.txt` were copied byte-for-byte from the supplied `Exactv2.2` directory on 2026-10-22. `original/SHA256.json` and the `original_hashes` CTest verify the import. Exact v1 is retained for provenance but is not an active comparison entry. The supplied README contains historical results, not results reproduced by this repository.

Call `exact_initialize()` once, serially, before using `exact_one` or `exact_batch` from `exact_adapter`. The boolean first argument selects the vectorized kernel. These wrappers validate inputs, map direction conventions and reject multi-revolution requests (`-6`). They do not extend the original kernels or repair inaccurate finite answers. Batched arbitrary 3D inputs are rotated into local transfer planes before the original shared-plane API is called, and velocities are rotated back. That work is included in comparison timing.

No license declaration was supplied with these Exact source files; this package does not invent one. No ivLam source, coefficient data, compiled binaries or external submodules are required or bundled.

See [validation](docs/VALIDATION.md) for measured results and limits. The repository is [Ossama-Abdelkhalik-Lab/eLam](https://github.com/Ossama-Abdelkhalik-Lab/eLam). It can be developed independently or checked out as the `solvers/elam` submodule of Lamberts-Comparison; both use the same Git history. Commit and push solver changes here first, then record the updated submodule commit in the comparison repository.
