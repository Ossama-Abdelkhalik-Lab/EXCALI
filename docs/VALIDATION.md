# Validation — 2026-10-22

This file records local evidence from Lamberts-Comparison on macOS with GNU Fortran 16.2 and MATLAB R2026b. It is not a claim of global superiority to ivLam or a reproduction of the supplied five-million-case benchmark.

## Completed checks

- Full native comparison suite: **24/24 CTests passed**.
- Independent standalone repository: **3/3 tests passed**, with OpenMP enabled and disabled.
- MATLAB R2026b: `build_solvers`, `run_tests`, `test_elam_comparison` and `test_comparison` all completed successfully, including report/plot export and replay. MATLAB required execution outside the filesystem sandbox to initialize its Qt runtime reliably.
- Original Exact files were also compared directly against the supplied directory, byte-for-byte. Russell's existing vendor-hash checks passed.

## Numerical coverage

- Preserved serial/vector originals: 517 independently known circular transfers in arbitrary 3D planes, both transfer directions, scalar/batch paths, and explicit unsupported multi-revolution/hyperbolic rejection. The import hash test checks all preserved files.
- Extended eLam: independent circular solutions through 10 revolutions, analytic parabolic/hyperbolic states, exact antipodal solutions with explicit transfer planes, invalid input and missing-plane statuses, and length/mu scales spanning `1e-150/1e-300` to `1e150/1e300`.
- 50,000 deterministic mixed requests compared against confluentcontour: elliptic, short/long hyperbolic, near-parabolic and both multi-revolution branches. This is cross-solver agreement, separate from independent truth checks.
- Near-antipodal and near-coincident/full-period cases checked by independent two-body propagation, including angular gaps down to approximately `1e-12` radians. Near-antipodal cases include both branches through 10 revolutions and extreme unit scaling.
- 32,400 multi-revolution feasibility decisions against a separately minimized time function, for 1–200 revolutions and relative distances `1e-1` through `1e-6` above/below the minimum. Feasible cases must return a success/merged status; failures cannot pass merely by avoiding status 5.
- MATLAB MEX integration: source/vendor hashes, MATLAB/native contour parity, dimensional scaling, initialization, input validation, reports, imported scalar and rotated-batch paths, and exact one/four-thread repeatability.

The MATLAB hard mixed fixture (`seed=173`) contains 3,000 generated cases plus nine appended stress cases. eLam passes all 3,000 generated cases and four stress cases: **3,004/3,009**, matching confluentcontour's pass count and failing case IDs. The same five stress cases remain visible: three exceed the endpoint velocity tolerance and two are singular without a transfer plane. Russell passes 3,003/3,009. Each preserved original passes 788/3,009; unsupported cases remain in the denominator. These totals do not establish the relative accuracy of the solvers within the originals' supported domain.

The comparison gate requires both endpoint velocities within `1e-10 + 1e-10*norm(truthVelocity)`. The near-antipodal propagation tests use a tighter normalized terminal-state tolerance of `1e-11`. The mixed cross-solver test uses regime-specific thresholds, including `2e-9` for strongly conditioned long-way hyperbolas.

## Timing interpretation

The numbers in this section and the next were measured before the zero-height rectangle update (see the end of this file). At that time the extension was slower than ivLam on the local mixed benchmark. Original and extended timings must be reported separately.

Original batched adapters include input validation, rotation into individual local transfer planes and rotation of output velocities back to the input frame. Serial and vector variants retain their supplied kernel-specific compiler flags. The common build uses `-O3 -fno-fast-math`; no `-march=native` speed claim is made. Initialization is outside timed regions. Thread counts are matched, and the parallel build is also checked at one thread.

The supplied `original/README.txt` contains historical timing and percentile claims from a different dataset and build. Those claims are retained as provenance, not adopted as results of this validation. No five-million-case corpus was supplied in the imported folder.

## Final native measurements

100,000 independent circular transfers; median of five repetitions; nanoseconds per complete solve. All five entries passed all 100,000 cases in each mode.

| Solver | Scalar, 1 thread | Batch, 1 thread | Batch, 4 threads |
|---|---:|---:|---:|
| confluentcontour | 236.3 | 134.4 | 36.6 |
| ivLam | 236.0 | 228.1 | 71.0 |
| Original Exact serial | 253.3 | 90.4 | 25.3 |
| Original Exact vector | 269.7 | 94.7 | 25.1 |
| Extended eLam | 1005.8 | 995.2 | 289.8 |

Thus the originals were about **2.4–2.5x faster in one-thread batches** and **2.8x faster in four-thread batches** on this circular fixture. Their scalar adapters were slower than ivLam. Maximum absolute velocity errors on this fixture were `2.68e-15` for ivLam, `3.16e-15`/`3.76e-15` for the originals and `4.27e-15` for extended eLam; this does not establish superior original accuracy.

The 22,000-case independent mixed benchmark (2,000 per class) produced zero failures for extended eLam, confluentcontour and ivLam, including 2,000 correctly classified infeasible requests. The originals each failed five high-eccentricity cases, with a worst finite relative endpoint error of **0.4079**. They also rejected unsupported hyperbolic/multi-revolution requests and about half the near-parabolic class. Ordinary short-/long-way elliptic classes passed, but the originals' worst long-way relative error (`1.05e-11`) exceeded ivLam's (`3.44e-13`) on this sample. The statement "more accurate than ivLam" is therefore not supported as a general claim.

Extended eLam took roughly **4.2–6.1x** ivLam's time on the feasible mixed transfer classes, and 7.3x for infeasibility classification. Its finite-rectangle method has broader coverage, but speed optimization remains future work rather than a result claimed here.

Raw results: [one-thread batch](circular_batch_1thread.csv), [four-thread batch](circular_batch_4threads.csv), [one-thread scalar](circular_scalar_1thread.csv), [mixed regimes](extended_mixed.csv). [Metadata](benchmark_metadata.json) records sample sizes, statistics, compiler, solver IDs and the comparison base revision. [Source hashes](validated_sources.json) identify the validated eLam implementation. Mixed timings use the minimum of five repeated block timings; the circular table uses medians. Final native timings were run after the task's builds, propagation checks and MATLAB rendering finished.

## Reproduction

Standalone build/test commands are in the main README. Full comparison tests are in Lamberts-Comparison:

```sh
cmake -S . -B build/elam-validation -DCMAKE_BUILD_TYPE=Release
cmake --build build/elam-validation -j
ctest --test-dir build/elam-validation --output-on-failure
```

```matlab
build_solvers;
addpath('tests');
run_tests;
test_elam_comparison;
test_comparison;
```

Native `benchmark_mixed` retains independent quad-precision orbital truth and per-regime failure counts for IDs 0, 1, 6, 7, 8. `benchmark_native_solvers` exercises complete scalar/batch gateway paths on circular truth. Rejected/failed rows do not justify speed claims. See the CSV files accompanying this report for the final local measurements.

## Zero-height rectangle update — 2026-10-22

The root step now uses the zero-height limit of the thin rectangle in a fixed one- or two-stage pipeline; see the README's Method section. Bracketed Newton/bisection remains only as a fallback for rare starter misses.

- Lamberts-Comparison native CTests: **24/24 passed** (including near-antipodal and near-coincident propagation, 32,400 minimum-time feasibility decisions, and the 50,000-case cross-solver agreement test).
- Standalone repository: **3/3 passed**, with OpenMP on and off.

`benchmark_mixed`, 2,000 cases per class, five repetitions, one thread, Apple silicon, GNU Fortran 16.2, `-O3`. Errors are maximum relative endpoint velocity errors against quad-precision Kepler truth. All solvers had zero failures.

| Class | ivLam ns | eLam ns | eLam / ivLam | ivLam max error | eLam max error |
|---|---:|---:|---:|---:|---:|
| circular_short | 200 | 257 | 1.28 | 3.6e-15 | 1.6e-15 |
| elliptic_short | 199 | 279 | 1.40 | 1.1e-13 | 4.4e-14 |
| elliptic_long | 192 | 248 | 1.29 | 3.4e-13 | 2.7e-14 |
| high_ecc | 188 | 320 | 1.70 | 1.6e-12 | 1.0e-14 |
| hyperbolic | 198 | 268 | 1.35 | 2.2e-14 | 6.4e-15 |
| near_parabolic | 178 | 258 | 1.45 | 1.5e-15 | 2.0e-15 |
| multirev_1 | 256 | 447 | 1.75 | 8.4e-13 | 1.2e-13 |
| multirev_2_3 | 258 | 425 | 1.64 | 4.9e-13 | 1.4e-13 |
| multirev_long | 257 | 419 | 1.63 | 1.1e-12 | 4.4e-13 |
| infeasible_m | 52 | 183 | 3.52 | — | — |
| hyperbolic_fast | 203 | 241 | 1.19 | 2.9e-15 | 8.2e-15 |

The previous eLam implementation took 4.5–5.3 times ivLam's time on the same benchmark (7.0 for infeasible requests).

On the Lamberts-Comparison 1,000,000-case elliptic dataset, the update was checked on 95,748 selected cases: 50,000 random cases, every transfer within 10° of 0° or 360°, and every transfer within 1° of 180°. The maximum endpoint error matched the previous implementation in every group: 3.1e-13 for the random cases, 1.6e-14 near 360°, 1.2e-12 below 10°, and 6.9e-11 near 180°. Near 180°, 21 cases exceed 1e-12 (20 before); this is the input-conditioning floor that all solvers share there. Fallbacks to the bracketed iteration occurred in 0–2% of cases per group.

## Speed update — 2026-10-23

The root extraction is unchanged: closed-form starter, one 8-node zero-height rectangle, and a 3-node rectangle only when the first stage's own error estimate is above rounding. The speed-up comes from the surrounding work:

- **Infeasible multi-revolution times** are rejected by the bound `T > m*2*pi*(s/2)^1.5` (from `a >= s/2`) before `T` is evaluated. This bound lies within 31% of the true minimum time.
- **Minimum-time search:** the large-`m` starter now uses `y(pi)=y0+2b` for `b>0` (it previously used the `b<0` value), which cuts its error from up to 0.31 to about `0.15/m`. Halley stops once its measured error bounds (`12*step^3` in `v`, `300*step^4` relative in `tmin`) decide feasibility. Half-angle functions come from `cos v` and then by angle addition, so the search makes no trig calls. Within 1% of the minimum time, the branch starter inverts the local cubic Taylor model.
- **Short-way elliptic starters** had the same `y(pi)` error in the pole asymptote. Fixing it reduces first-pass rectangle misses on short arcs from 37% to about 12%.
- **Rectangle kernel:** one fused loop per coordinate over compile-time tables of mirror nodes, constant multiplies instead of divisions, two partial sums, and one reciprocal for the ratio and error estimate. An 8-node stage now costs about 35 ns (66 before) in the half-anomaly coordinate and 45 ns (120 before) in the parabolic series coordinate.
- **Library calls:** starters use a bit-pattern cube root with one Halley correction (relative error 2e-5; the contour sets the accuracy). Branch-point distances use a fast logarithm. Powers of two and exponents come from bit patterns (gfortran otherwise calls `__powidf2` and `frexp`). The complex `**2` in the X-plane distance is replaced by real arithmetic. `sinh(wb)` uses an identity, and the centre's half-angle functions are carried between stages and into the velocities.

All 24 Lamberts-Comparison CTests pass. The standalone repository passes 3/3 with OpenMP on and off, and also with `ELAM_PGO=USE`.

`benchmark_mixed`, 2,000 cases per class, minimum of five block timings in each of five interleaved runs, one thread, Apple M1 Pro, GNU Fortran 16.2, `-O3 -fno-fast-math` (no PGO). All solvers had zero failures. "Before" is eLam at commit 50b9585.

| Class | ivLam ns | eLam before | eLam now | now / ivLam | ivLam max error | eLam max error |
|---|---:|---:|---:|---:|---:|---:|
| circular_short | 198 | 253 | 158 | 0.80 | 3.6e-15 | 2.2e-15 |
| elliptic_short | 197 | 278 | 164 | 0.83 | 1.1e-13 | 4.4e-14 |
| elliptic_long | 186 | 242 | 169 | 0.91 | 3.4e-13 | 2.7e-14 |
| high_ecc | 185 | 318 | 185 | 1.00 | 1.6e-12 | 1.0e-14 |
| hyperbolic | 197 | 267 | 178 | 0.90 | 2.2e-14 | 6.4e-15 |
| near_parabolic | 176 | 253 | 158 | 0.90 | 1.5e-15 | 2.0e-15 |
| multirev_1 | 253 | 441 | 279 | 1.10 | 8.4e-13 | 5.8e-13 |
| multirev_2_3 | 253 | 424 | 248 | 0.98 | 4.9e-13 | 2.2e-13 |
| multirev_long | 254 | 411 | 247 | 0.97 | 1.1e-12 | 3.9e-13 |
| infeasible_m | 51 | 180 | 21 | 0.40 | — | — |
| hyperbolic_fast | 202 | 236 | 154 | 0.76 | 2.9e-15 | 6.8e-15 |
| **overall** | | | | **0.91** (1.53 before) | | |

Average work per solve, from instrumented builds:

| Class | 8-node stages | 3-node stage | Fallback |
|---|---:|---:|---:|
| circular_short | 1.13 | 2.2% | 0.20% |
| elliptic_short | 1.10 | 14% | 0.10% |
| elliptic_long | 1.02 | 14% | 0.25% |
| high_ecc | 1.16 | 17% | 0.70% |
| hyperbolic | 1.11 | 14% | 0.40% |
| near_parabolic | 1.00 | 0.5% | 0 |
| multirev_1 | 1.06 | 11% | 0.40% |
| multirev_2_3 | 1.04 | 6% | 0.15% |
| multirev_long | 1.05 | 5% | 0.10% |
| hyperbolic_fast | 1.00 | 0.2% | 0 |

**Accuracy:**
- **95,748 selected cases** from the 1,000,000-case elliptic dataset: every case is within 4× of the previous implementation's error (or below 1e-14). The group maxima are unchanged: 3.1e-13 for the random cases, 1.6e-14 near 360°, 1.2e-12 below 10° and 6.9e-11 near 180°.
- **Multi-revolution, 20,000 cases per class** under three seeds, compared case by case: p99.9 errors match the previous implementation within ±30%. Cases that got worse and cases that got better occur in similar numbers. All of them lie within about 1e-5 of the minimum time, where the root is poorly conditioned, and the two versions agree to about 1e-13. Per-seed worst cases moved both ways: 1.3e-9 → 1.0e-9 and 1.5e-10 → 1.6e-11 on two seeds, and 3.3e-12 → 1.6e-11 on the third.

**Profile-guided build:** `ELAM_PGO=USE` makes eLam a further 3–5% faster in every class. The table above does not use it, because the comparison builds every solver with the same flags.

**Measured on this machine and not adopted.** These are platform-dependent trade-offs. The adopted code has no platform-specific paths, so none of these choices is built into it:
- **Mixed precision:** a single-precision stage 1 saved 4–11 ns per 8-node stage. It then needs a double-precision stage 2 on every solve, instead of on about 10% of solves, so here it was a net loss. With wider single-precision vectors (for example AVX-512) the balance could change.
- **Derivative-based (Turán-type) quadrature:** per node, the kernel costs about 2 ns here (5 nodes against 10 differ by about 10 ns). Using `T''` to drop three nodes would save under 6 ns, before its extra work.
- **Inline polynomial sin/cos:** accurate to 1.4 ulp, but no faster than the library call in the real solve. The library pair is kept everywhere. A platform-specific entry such as Apple's `__sincos_stret` gained about 1% and was rejected as nonportable.
- **`-mcpu=native`, `-funroll-loops`, `-fno-math-errno`:** no measurable effect here.
- **Batching problems in lockstep:** not attempted. On this machine's two-double vectors, the remaining time is a dependency chain within each solve rather than arithmetic throughput. Wider x86 vectors could make it worthwhile.

### Portability

The solver uses only standard Fortran 2008 and the C99 functions `cbrt` and `log1p`. The fast cube root, logarithm, powers of two and exponents read the IEEE-754 binary64 bit pattern. A compile-time check (`ieee64`) switches them to the standard intrinsics on any other real format. A build with that switch forced off passed the same accuracy checks.

Also checked in Docker (`gcc:14` image, GNU Fortran 14.4, glibc):

- **Linux aarch64 and Linux x86-64:** the standalone tests pass, 3/3 with OpenMP on and 3/3 with it off. The x86-64 run is emulated, so it gives no timing information.
- **95,748 selected cases on Linux aarch64:** the errors are identical to macOS.
- **95,748 selected cases on Linux x86-64:** no case is worse than 4× the previous eLam on that platform, and every group maximum is unchanged from it. In the near-180° group, both versions differ from their macOS results in 47 cases, with a maximum of 1.0e-10 against 6.9e-11. That group sits at the input-conditioning floor, where x86-64's different rounding (no fused multiply-add by default) and different libm show.

All the timings in this file are from one Apple M1 Pro. Native x86-64 timings have not been measured.
