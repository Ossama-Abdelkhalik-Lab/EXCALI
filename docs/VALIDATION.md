# Validation — 2026-10-22

This file records local evidence from Lamberts-Comparison on macOS with GNU Fortran 16.2 and MATLAB R2026b. It is not a claim of global superiority to ivLam or a reproduction of the supplied five-million-case benchmark.

## Completed checks (2026-10-22)

- Full native comparison suite: **24/24 CTests passed**.
- Independent standalone repository: **3/3 tests passed**, with OpenMP enabled and disabled.
- MATLAB R2026b: `build_solvers`, `run_tests`, `test_excali_comparison` and `test_comparison` all completed successfully, including report/plot export and replay. MATLAB required execution outside the filesystem sandbox to initialize its Qt runtime reliably.
- Original Exact files were also compared directly against the supplied directory, byte-for-byte. Russell's existing vendor-hash checks passed.

## Numerical coverage

- Preserved serial/vector originals: 517 independently known circular transfers in arbitrary 3D planes, both transfer directions, scalar/batch paths, and explicit unsupported multi-revolution/hyperbolic rejection. The import hash test checks all preserved files.
- Extended EXCALI: independent circular solutions through 10 revolutions, analytic parabolic/hyperbolic states, exact antipodal solutions with explicit transfer planes, invalid input and missing-plane statuses, and length/mu scales spanning `1e-150/1e-300` to `1e150/1e300`.
- 50,000 deterministic mixed requests compared against confluentcontour: elliptic, short/long hyperbolic, near-parabolic and both multi-revolution branches. This is cross-solver agreement, separate from independent truth checks.
- Near-antipodal and near-coincident/full-period cases checked by independent two-body propagation, including angular gaps down to approximately `1e-12` radians. Near-antipodal cases include both branches through 10 revolutions and extreme unit scaling.
- 32,400 multi-revolution feasibility decisions against a separately minimized time function, for 1–200 revolutions and relative distances `1e-1` through `1e-6` above/below the minimum. Feasible cases must return a success/merged status; failures cannot pass merely by avoiding status 5.
- MATLAB MEX integration: source/vendor hashes, MATLAB/native contour parity, dimensional scaling, initialization, input validation, reports, imported scalar and rotated-batch paths, and exact one/four-thread repeatability.

The MATLAB hard mixed fixture (`seed=173`) contains 3,000 generated cases plus nine appended stress cases. EXCALI passes all 3,000 generated cases and four stress cases: **3,004/3,009**, matching confluentcontour's pass count and failing case IDs. The same five stress cases remain visible: three exceed the endpoint velocity tolerance and two are singular without a transfer plane. Russell passes 3,003/3,009. Each preserved original passes 788/3,009; unsupported cases remain in the denominator. These totals do not establish the relative accuracy of the solvers within the originals' supported domain.

The comparison gate requires both endpoint velocities within `1e-10 + 1e-10*norm(truthVelocity)`. The near-antipodal propagation tests use a tighter normalized terminal-state tolerance of `1e-11`. The mixed cross-solver test uses regime-specific thresholds, including `2e-9` for strongly conditioned long-way hyperbolas.

## Timing interpretation

Original and extended timings must be reported separately. The original batched adapters include input validation, rotation into individual local transfer planes and rotation of output velocities back to the input frame. Serial and vector variants retain their supplied kernel-specific compiler flags. The common build uses `-O3 -fno-fast-math`; no `-march=native` speed claim is made. Initialization is outside timed regions. Thread counts are matched, and the parallel build is also checked at one thread.

The supplied `original/README.txt` contains historical timing and percentile claims from a different dataset and build. Those claims are retained as provenance, not adopted as results of this validation. No five-million-case corpus was supplied in the imported folder.

On the 22,000-case mixed benchmark (2,000 per class) the preserved originals each failed five high-eccentricity cases, with a worst finite relative endpoint error of 0.4079, and rejected the unsupported hyperbolic and multi-revolution requests and about half the near-parabolic class. Their worst long-way elliptic error (`1.05e-11`) exceeded ivLam's (`3.44e-13`), so "more accurate than ivLam" is not supported as a general claim. Timings of the first finite-rectangle EXCALI, before the zero-height update below, are superseded and no longer recorded here.

## Reproduction

Standalone build/test commands are in the main README. Full comparison tests are in Lamberts-Comparison:

```sh
cmake -S . -B build/excali-validation -DCMAKE_BUILD_TYPE=Release
cmake --build build/excali-validation -j
ctest --test-dir build/excali-validation --output-on-failure
```

```matlab
build_solvers;
addpath('tests');
run_tests;
test_excali_comparison;
test_comparison;
```

Native `benchmark_mixed` retains independent quad-precision orbital truth and per-regime failure counts for confluentcontour, Russell, both preserved originals and EXCALI (comparison IDs 0–4 since 2026-10-24). `benchmark_native_solvers` exercises complete scalar/batch gateway paths on circular truth. Rejected/failed rows do not justify speed claims.

## Zero-height rectangle update — 2026-10-22

The root step now uses the zero-height limit of the thin rectangle in a fixed one- or two-stage pipeline; see the README's Method section. Bracketed Newton/bisection remains only as a fallback for rare starter misses.

- Lamberts-Comparison native CTests: **24/24 passed** (including near-antipodal and near-coincident propagation, 32,400 minimum-time feasibility decisions, and the 50,000-case cross-solver agreement test).
- Standalone repository: **3/3 passed**, with OpenMP on and off.

`benchmark_mixed`, 2,000 cases per class, five repetitions, one thread, Apple silicon, GNU Fortran 16.2, `-O3`. Errors are maximum relative endpoint velocity errors against quad-precision Kepler truth. All solvers had zero failures.

| Class | ivLam ns | EXCALI ns | EXCALI / ivLam | ivLam max error | EXCALI max error |
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

The previous EXCALI implementation took 4.5–5.3 times ivLam's time on the same benchmark (7.0 for infeasible requests).

On the Lamberts-Comparison 1,000,000-case elliptic dataset, the update was checked on 95,748 selected cases: 50,000 random cases, every transfer within 10° of 0° or 360°, and every transfer within 1° of 180°. The maximum endpoint error matched the previous implementation in every group: 3.1e-13 for the random cases, 1.6e-14 near 360°, 1.2e-12 below 10°, and 6.9e-11 near 180°. Near 180°, 21 cases exceed 1e-12 (20 before); this is the input-conditioning floor that all solvers share there. Fallbacks to the bracketed iteration occurred in 0–2% of cases per group.

## Speed update — 2026-10-23

The root extraction is unchanged: closed-form starter, one 8-node zero-height rectangle, and a 3-node rectangle only when the first stage's own error estimate is above rounding. The speed-up comes from the surrounding work:

- **Infeasible multi-revolution times** are rejected by the bound `T > m*2*pi*(s/2)^1.5` (from `a >= s/2`) before `T` is evaluated. This bound lies within 31% of the true minimum time.
- **Minimum-time search:** the large-`m` starter now uses `y(pi)=y0+2b` for `b>0` (it previously used the `b<0` value), which cuts its error from up to 0.31 to about `0.15/m`. Halley stops once its measured error bounds (`12*step^3` in `v`, `300*step^4` relative in `tmin`) decide feasibility. Half-angle functions come from `cos v` and then by angle addition, so the search makes no trig calls. Within 1% of the minimum time, the branch starter inverts the local cubic Taylor model.
- **Short-way elliptic starters** had the same `y(pi)` error in the pole asymptote. Fixing it reduces first-pass rectangle misses on short arcs from 37% to about 12%.
- **Rectangle kernel:** one fused loop per coordinate over compile-time tables of mirror nodes, constant multiplies instead of divisions, two partial sums, and one reciprocal for the ratio and error estimate. An 8-node stage now costs about 35 ns (66 before) in the half-anomaly coordinate and 45 ns (120 before) in the parabolic series coordinate.
- **Library calls:** starters use a bit-pattern cube root with one Halley correction (relative error 2e-5; the contour sets the accuracy). Branch-point distances use a fast logarithm. Powers of two and exponents come from bit patterns (gfortran otherwise calls `__powidf2` and `frexp`). The complex `**2` in the X-plane distance is replaced by real arithmetic. `sinh(wb)` uses an identity, and the centre's half-angle functions are carried between stages and into the velocities.

All 24 Lamberts-Comparison CTests pass. The standalone repository passes 3/3 with OpenMP on and off, and also with `EXCALI_PGO=USE`.

`benchmark_mixed`, 2,000 cases per class, minimum of five block timings in each of five interleaved runs, one thread, Apple M1 Pro, GNU Fortran 16.2, `-O3 -fno-fast-math` (no PGO). All solvers had zero failures. "Before" is EXCALI at commit 50b9585.

| Class | ivLam ns | EXCALI before | EXCALI now | now / ivLam | ivLam max error | EXCALI max error |
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

**Profile-guided build:** `EXCALI_PGO=USE` makes EXCALI a further 3–5% faster in every class. The table above does not use it, because the comparison builds every solver with the same flags.

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
- **95,748 selected cases on Linux x86-64:** no case is worse than 4× the previous EXCALI on that platform, and every group maximum is unchanged from it. In the near-180° group, both versions differ from their macOS results in 47 cases, with a maximum of 1.0e-10 against 6.9e-11. That group sits at the input-conditioning floor, where x86-64's different rounding (no fused multiply-add by default) and different libm show.

Every timing above is from one Apple M1 Pro.

### Native Linux x86-64 (GitHub Actions)

The comparison repository's `native` workflow runs `benchmark_mixed` (2,000 cases per class, five repetitions) on an `ubuntu-latest` runner: x86-64, GNU Fortran 13.3, the same `-O3 -fno-fast-math` flags. All solvers had zero failures. Shared runners are noisy, so the ratios matter more than the absolute times.

| Class | ivLam ns | EXCALI ns | EXCALI / ivLam | ivLam max error | EXCALI max error |
|---|---:|---:|---:|---:|---:|
| circular_short | 264 | 287 | 1.09 | 3.5e-15 | 1.8e-15 |
| elliptic_short | 264 | 321 | 1.22 | 1.5e-13 | 4.3e-14 |
| elliptic_long | 255 | 295 | 1.16 | 2.9e-13 | 4.1e-14 |
| high_ecc | 259 | 346 | 1.34 | 2.1e-12 | 1.9e-14 |
| hyperbolic | 271 | 379 | 1.40 | 3.1e-14 | 7.8e-15 |
| near_parabolic | 236 | 328 | 1.39 | 1.6e-15 | 4.9e-15 |
| multirev_1 | 331 | 434 | 1.31 | 8.1e-13 | 6.7e-14 |
| multirev_2_3 | 330 | 391 | 1.19 | 5.0e-13 | 2.0e-13 |
| multirev_long | 332 | 391 | 1.18 | 5.4e-13 | 4.1e-13 |
| infeasible_m | 93 | 57 | 0.61 | — | — |
| hyperbolic_fast | 280 | 301 | 1.07 | 4.7e-15 | 1.0e-14 |

A follow-up profile on one x86-64 runner compared the previous and current EXCALI in the same job. Computing sinh/cosh from one `exp` (or a Taylor series), and removing redundant divisions and square roots, changed the ratio to ivLam as follows:
- hyperbolic: 1.41 → 1.31;
- hyperbolic_fast: 1.07 → 1.00;
- near_parabolic: 1.39 → 1.35;
- all classes combined: 1.21 → 1.18.

The baseline x86-64 target used by plain `-O3` has no fused multiply-add, which EXCALI's polynomial evaluations rely on more than ivLam does.

So EXCALI's speed relative to ivLam depends on the platform. It is 0.76–1.10× on Apple silicon with macOS's libm, and 1.07–1.40× on this x86-64 runner. Its accuracy advantage holds on both. The same pattern appeared on Linux aarch64 in Docker on the M1: glibc's `sinh`/`cosh`, `cbrt` and `asinh` took about twice Apple's time, and the rectangle kernel ran more slowly. Those Docker timings ran inside a VM, though, so they are only indicative.

## Contour-only safeguards — 2026-10-24

The two derivative-based iterations left in the root solve were replaced with zero-height rectangles:

- **Starter misses** used to fall back to bracketed Newton/bisection, whose first step was a finite complex rectangle. They now continue as a *bracketed contour*. Every estimate is a zero-height rectangle ratio centred on the current point. `T` is monotone on the branch, so two pieces of information shrink the starter's bracket: the sign of `T - t` at the centre, and the rectangle's own sign test (a sign change confines the root to the rectangle; none moves the bracket edge past it). An estimate outside the bracket is replaced by its midpoint. A rectangle that leaves the domain (for example past the hyperbolic branch point) returns a non-finite ratio and only halves `L`.
- **Very short hyperbolas** (`t < 0.1 t_p`) solved `T(Y) = t` in `Y = sqrt(y)` by Newton with a complex-step derivative. They now use the same bracketed contour in `Y`, with the closed-form `dT/dY = 3Y^2 F + a - 2Y^4 F'(w)/(b sinh w)`. The anomaly coordinate `w` cannot replace `Y` there: in a stress test of 20,000 hyperbolas with `e` from 30 to 10,000, 72% of the cases failed in `w` because the root sits within rounding of the branch point.

The finite complex rectangle and the complex time function are gone. EXCALI no longer uses complex arithmetic anywhere. The only derivative-based iteration left is the Halley search for the multi-revolution minimum time (one or two steps).

Removing the safeguard altogether is not possible without failures. With it disabled, 0.02–0.77% of each class of `benchmark_mixed` return status `-4` (20,000 cases per class, eight seeds). In most of these the starter is far from the root and the reach-limited interval never contains it (for example a long-way hyperbola whose `T - t` stays near +0.31 across the whole interval), or the root sits beside a pole of `T`.

**Checks:**

- `benchmark_mixed` with 20,000 cases per class under eight seeds (1.76 million solves): zero failures. Each class maximum matches the previous Newton fallback; the rescued cases take 3–6 rectangles each.
- Short-hyperbola stress set, five seeds of 20,000 (75,000 cases on the `Y` path): zero failures, maxima equal to the Newton version (2.4e-14 to 3.2e-14).
- The full 1,000,000-case elliptic dataset, case by case against the previous build: 999,930 results are bit-identical, no case differs by more than 4×, there are no failures, and the maxima are unchanged (1.2e-13 for 10°–179°, 3.4e-13 below 10°, 8.8e-11 above 179°).
- Lamberts-Comparison CTests: 24/24 passed. Standalone repository: 3/3 with OpenMP on and off. MATLAB `run_tests` and `test_excali_comparison`: passed; the hard mixed fixture remains 3,004/3,009.
- `benchmark_mixed` speed is unchanged within noise (0.87 times ivLam overall on the M1 Pro).

**Rounding sensitivity (found by the Linux x86-64 CI).** Plain `-O3` on x86-64 does not fuse multiply-adds, and a Mac build with `-ffp-contract=off` reproduces its rounding. Under that rounding the first version of the bracketed contour failed 653 near-parabolic and one near-minimum multi-revolution case. Three rules fix it:
- near the parabola its rectangles use the series coordinate `X = ±x²`, as the main stages do, because the closed form in `x` cancels;
- a rectangle without a sign change moves the bracket only when its own ratio also puts the root outside, since near the root both signals are rounding noise;
- rectangles stay inside the bracket as it was before each step, so they never reach the other branch beyond `t_min`.

With these, both roundings match the previous Newton fallback in every `benchmark_mixed` class over eight seeds (same failures and maxima; the one no-FMA `multirev_2_3` case at 1.2e-10 is shared with the old code), and 999,937 of the 1,000,000 dataset results are bit-identical. The same runs show that without FMA the main stages miss far more often: 308–3,472 of 160,000 cases per class instead of 0–1,133. More solves therefore reach the bracketed contour on x86-64, which is a likely contributor to its slower timings.

## Starter tables, one-rectangle path and batch entry — 2026-10-27

Compiled-in starter tables (`src/excali_tables.inc`, generated by `tools/gen_tables.py`) now give the root for zero-revolution ellipses and hyperbolas, and for multiple revolutions the branch root, the minimum time and its location. A second-order series reversion covers the near-parabolic band outside the tables. From these starters one 8-node rectangle, with half-width at most 1/16 of a cheap lower bound on the singularity distance, is final; any rectangle it does not certify (sign change, estimate inside, error estimate at rounding) continues on the previous general path. Multi-revolution feasibility comes from the table's `t_min` outside a 1e-6 margin; inside it, the Halley search runs as before. `excali_batch` runs the common path stage by stage over blocks of 32 cases. Other changes: an inline half-angle sine/cosine (the `sin`/`cos` pair compiled to a `cexp` call), and a cheaper polynomial logarithm for branch-point distances. Tried and not kept: putting the node formulas over common denominators (one division per node instead of three; slower), and `-ffast-math` (5% faster, but two high-eccentricity failures and two crashing tests).

**Checks:**
- Standalone repository: 3/3 with OpenMP on and off, and with `EXCALI_PGO=USE`.
- Lamberts-Comparison: 9/9 CTests; MATLAB `build_solvers`, `test_excali_comparison` (3,004/3,009 on the hard mixed fixture, unchanged, including 1/4-thread parity of the batch path), `test_plots` and `run_tests`.
- 400,000 extreme random requests (radius ratios to 1e3, transfer angles within 1e-8 of 0, pi and 2 pi, times from 1e-6 to 1e6 times the transfer scale, up to 20 revolutions on both branches) against the previous implementation: identical statuses for every case (394,232 solved, 5,755 infeasible, 13 unsolved by both), scalar and batch. All solutions agree within 1e-11 except one at 2.6e-10, a radius ratio of 1e-3 with a 4e-7 rad transfer angle, where independent propagation puts both solutions within 3.4e-11 of the target.
- 1,000,009-case hard mixed dataset (seed 4242): identical pass counts and statuses to the previous implementation in every category.
- Fast-path coverage (20,000 cases per class): the one-rectangle path finishes 99.4–100% of the cases the starters cover; 94.5% of the `hyperbolic` class is covered (the rest are very short hyperbolas, `t < 0.1 t_p`, which keep the bracketed contour in `sqrt(y)`), and 99.8% of multi-revolution cases.

`benchmark_mixed`, 2,000 cases per class, minimum of seven block timings, one thread, Apple M1 Pro, GNU Fortran 16.2, `-O3 -fno-fast-math`. The machine was under background load; ratios within a run are reliable, absolute times to about ±5%. "Before" is EXCALI at commit 4f759c8 in a run with ivLam at 198 ns overall.

| Class | ivLam ns | confluentcontour ns | Exact serial ns | EXCALI before | EXCALI now | now / ivLam |
|---|---:|---:|---:|---:|---:|---:|
| circular_short | 199 | 187 | 213 | 160 | 150 | 0.75 |
| elliptic_short | 198 | 198 | 220 | 169 | 152 | 0.77 |
| elliptic_long | 187 | 184 | 219 | 173 | 153 | 0.82 |
| high_ecc | 188 | 194 | 227 (5 fail) | 183 | 151 | 0.80 |
| hyperbolic | 198 | 217 | rejects | 177 | 164 | 0.83 |
| near_parabolic | 176 | 138 | 200 (1025 fail) | 159 | 135 | 0.76 |
| multirev_1 | 256 | 203 | rejects | 281 | 167 | 0.65 |
| multirev_2_3 | 257 | 200 | rejects | 252 | 167 | 0.65 |
| multirev_long | 257 | 234 | rejects | 253 | 180 | 0.70 |
| infeasible_m | 51 | 46 | 16 | 26 | 24 | 0.47 |
| hyperbolic_fast | 201 | 184 | rejects | 152 | 157 | 0.78 |
| **overall** | 197 | 180 | | 181 (0.91) | **146** | **0.74** |

The batch entry, same classes and cases (native harness, interleaved with the single-solve loop of the previous implementation, ns per case): circular_short 90, elliptic_short 94, elliptic_long 101, high_ecc 105, hyperbolic 117, near_parabolic 105, multirev_1 116, multirev_2_3 117, multirev_long 137, infeasible_m 23, hyperbolic_fast 110; 101 overall against 182.

On the 1,000,000-case zero-revolution elliptic dataset of the MATLAB comparison (native batch routines, one thread, interleaved): EXCALI `excali_batch` 118 ns per case, Exact v2.2 serial batch 133, Exact v2.2 vector batch 120. EXCALI fails none; each Exact build fails 3 at the MATLAB tolerance. On the hard mixed dataset Exact's batches are faster per case only because they reject the 74% of cases outside their domain.

Accuracy on 20,000 cases per class against quad-precision truth (relative endpoint velocity error):

| Class | EXCALI p99 | confluentcontour p99 | ivLam p99 | EXCALI p99.9 | confluentcontour p99.9 | ivLam p99.9 |
|---|---:|---:|---:|---:|---:|---:|
| circular_short | 1.2e-15 | 1.0e-15 | 2.0e-15 | 1.6e-15 | 1.4e-15 | 3.1e-15 |
| elliptic_short | 2.2e-15 | 2.4e-15 | 8.1e-15 | 8.0e-15 | 8.2e-15 | 6.0e-14 |
| elliptic_long | 1.8e-15 | 2.2e-15 | 5.1e-14 | 5.0e-15 | 5.5e-15 | 2.6e-13 |
| high_ecc | 4.7e-15 | 5.1e-15 | 7.5e-14 | 9.2e-15 | 9.8e-15 | 6.8e-13 |
| hyperbolic | 1.4e-15 | 1.6e-15 | 4.5e-15 | 4.6e-15 | 3.8e-15 | 2.1e-14 |
| near_parabolic | 1.2e-15 | 1.2e-15 | 7.9e-16 | 1.9e-15 | 2.0e-15 | 1.1e-15 |
| multirev_1 | 2.0e-14 | 2.0e-14 | 7.2e-14 | 1.8e-13 | 1.6e-13 | 3.5e-13 |
| multirev_2_3 | 1.9e-14 | 2.1e-14 | 6.8e-14 | 1.7e-13 | 2.2e-13 | 3.5e-13 |
| multirev_long | 1.9e-14 | 2.0e-14 | 7.2e-14 | 1.8e-13 | 1.9e-13 | 3.1e-13 |
| hyperbolic_fast | 4.0e-15 | 4.7e-15 | 1.2e-15 | 5.6e-15 | 7.5e-15 | 1.9e-15 |

ivLam fails 4 of the 60,000 multi-revolution cases; EXCALI and confluentcontour fail none. EXCALI's worst multi-revolution case (1.75e-11, unchanged from the previous implementation) lies 7.7e-11 above the minimum time, where the two branch roots are about 7e-6 apart; confluentcontour reaches 8e-13 there and ivLam fails.

**Hard stress cases.** Of the nine appended stress cases, `circular`, `high_eccentricity`, `near_parabolic` and `near_half_short` pass. The exact solution of the rounded double inputs, computed in quad precision, already misses the generating orbit by 1.6e-9 for `near_collinear` and 3e-10 for `near_half_long` (transfer angle within 1e-7 of 0 or pi: the orbit plane is known only to about 1e-9 from `r1 x r2`); `near_full_orbit` is similarly conditioned. `singular_half` and `singular_coincident` have no defined transfer plane. These five cannot meet the 1e-10 velocity tolerance from double inputs.

**Remaining iteration and fallbacks.** Very short hyperbolas (`t < 0.1 t_p`, 40% of the hard dataset's hyperbolas) still use the bracketed contour in `sqrt(y)` (1.2 rectangles per solve on average); multi-revolution times within 1e-6 of `t_min` still use the Halley search; and the general path still handles the 0–0.6% of covered cases whose single rectangle is rejected. Single solves of strong hyperbolas (`hyperbolic_fast`) are about 5% slower than before, because the table lookup costs more than the closed-form starter that already served them well; the batch entry is 27% faster there.

## Accuracy at the rounding floor, one code path — 2026-10-27 (later)

- **Velocities:** near the parabola and near the short-way branch point (both with cos(theta/2) > 0.3), Lagrange f,g on the exact inputs in double-double arithmetic (FMA products, compensated sums), correctly rounded for the computed `y`. On short-way hyperbolas near the branch point, `y` comes from the time equation, `sqrt(y) = t/(y*F + a)`, which damps the cancellation in `y0 - 2b*sinh(w/2)**2` (only where `a` dominates `y*F`; near 180 degrees it does not), and the half-angle functions follow from `y`. When the radii are within a factor 2, `|r1|-|r2|` comes from `(r1-r2).(r1+r2)/(|r1|+|r2|)`.
- **One code path:** `excali_one` is one case of the batch stages (`solve_chunk`), so geometry, table starters, the one-rectangle path and velocities each exist once; the general path is `general_root`. The zero-revolution elliptic table is accurate enough for a 3-node rectangle (half-width 1/256 of the bound).

Accuracy on 20,000 cases per class against quad-precision truth (relative endpoint velocity error, p99 / p99.9 / max):

| Class | EXCALI | confluentcontour | ivLam |
|---|---|---|---|
| circular_short | 1.2e-15 / 1.7e-15 / 2.4e-15 | 1.0e-15 / 1.4e-15 / 1.9e-15 | 2.0e-15 / 3.1e-15 / 4.2e-15 |
| elliptic_short | 2.1e-15 / 7.0e-15 / 4.8e-14 | 2.4e-15 / 8.2e-15 / 4.8e-14 | 8.1e-15 / 6.0e-14 / 4.5e-13 |
| elliptic_long | 1.8e-15 / 4.7e-15 / 3.1e-14 | 2.2e-15 / 5.5e-15 / 3.1e-14 | 5.1e-14 / 2.6e-13 / 1.4e-12 |
| high_ecc | 4.8e-15 / 8.4e-15 / 2.9e-14 | 5.1e-15 / 9.8e-15 / 2.9e-14 | 7.5e-14 / 6.8e-13 / 2.8e-12 |
| hyperbolic | 1.2e-15 / 5.0e-15 / 3.9e-14 | 1.6e-15 / 3.8e-15 / 4.4e-14 | 4.5e-15 / 2.1e-14 / 1.0e-13 |
| near_parabolic | 6.8e-16 / 9.8e-16 / 1.5e-15 | 1.2e-15 / 2.0e-15 / 2.9e-15 | 7.9e-16 / 1.1e-15 / 1.6e-15 |
| multirev_1 | 1.9e-14 / 1.8e-13 / 1.8e-12 | 2.0e-14 / 1.6e-13 / 3.1e-12 | 7.2e-14 / 3.5e-13 / 2.6e-12 (1 fail) |
| multirev_2_3 | 2.0e-14 / 1.7e-13 / 1.6e-12 | 2.1e-14 / 2.2e-13 / 1.2e-12 | 6.8e-14 / 3.5e-13 / 3.7e-12 (1 fail) |
| multirev_long | 1.9e-14 / 1.6e-13 / 1.8e-11 | 2.0e-14 / 1.9e-13 / 2.2e-12 | 7.2e-14 / 3.1e-13 / 1.3e-12 (2 fail) |
| hyperbolic_fast | 1.0e-15 / 1.7e-15 / 2.8e-15 | 4.7e-15 / 7.5e-15 / 1.4e-14 | 1.2e-15 / 1.9e-15 / 3.0e-15 |

EXCALI is below ivLam at p99 and p99.9 in every class and at the maximum in every class except multirev_long (one case 7.7e-11 above the minimum time, where ivLam fails two cases). Medians are at the rounding floor for all three (2-5e-16; ivLam's is lower by up to 1.3e-16 in five classes). The quad-precision exact solution of the rounded inputs, rounded to double, sets the floor at about 6e-16 (p99) for near_parabolic and 9e-16 for hyperbolic_fast.

**Speed (batch routines, one thread, interleaved, Apple M1 Pro under background load):** on the 1,000,000-case zero-revolution elliptic dataset EXCALI 114 ns per case, Exact v2.2 vector 120, serial 133 (EXCALI 0 failures, Exact 3). On the zero-revolution ellipses of the hard mixed dataset (the cases Exact attempts, 66,667): EXCALI 115 ns, Exact vector 116, serial 131 (failures 2 versus 34). Over the whole hard dataset Exact's batches average 143-165 ns against EXCALI's 195 only because they reject 737,000 of the 1,000,009 cases; EXCALI solves all but 7. Single solves through `excali_one` cost about 8% more than before the code paths were merged (the price of one code path); they remain 0.88 times ivLam's time overall in `benchmark_mixed` (0.73-0.92 per class, except 1.02 for hyperbolic and 1.07 for hyperbolic_fast, where the double-double velocities add their cost).

**Stress:** 400,000 extreme random requests give identical statuses to the implementation at the start of the day (commit 36228ca). The 1,172 solutions that differ by more than 1e-11 are all within 1e-4 rad of a 180-degree transfer, where independent propagation shows neither version consistently closer (both miss by about 1e-9).

**Very short hyperbolas (same day, later).** A rectangle in `xi = wb - w`, the distance from the short-way branch point, with `y = 2b*sinh(wb - xi/2)*sinh(xi/2)` (full relative precision as `y -> 0`) and a table of `xi/(wb*u**2)` for `u = t/tp < 0.1`, replaces the bracketed `sqrt(y)` contour for about 90% of those cases. The hyperbolic tables moved to a grid uniform in `sqrt|tau|` (worst error 0.4% of the singularity bound, against 0.9%), which extends them to `|tau| = 4e-4`, transfers within about 0.1 degree of 180 degrees. Exact parabolas no longer take the double-double velocity path. On the 1,000,009-case hard mixed dataset the batch now takes about 152 ns per case (195 before these two changes); Exact v2.2's vector batch takes about 145 ns and its serial batch 166 ns while rejecting 74% of the cases (its vector kernel runs the whole pipeline on every case, so a rejected case costs about as much as a solved one). On the 1,000,000-case zero-revolution elliptic dataset EXCALI takes 115-117 ns per case against Exact's 121-124 (vector) and 134-137 (serial).

Maximum error over the cases both solvers solve (20,000 per class; ivLam fails 4 multi-revolution cases that EXCALI solves): EXCALI / ivLam circular_short 2.4e-15 / 4.2e-15, elliptic_short 4.8e-14 / 4.5e-13, elliptic_long 3.1e-14 / 1.4e-12, high_ecc 2.9e-14 / 2.8e-12, hyperbolic 3.9e-14 / 1.0e-13, near_parabolic 1.5e-15 / 1.6e-15, multirev_1 1.8e-12 / 2.6e-12, multirev_2_3 1.6e-12 / 3.7e-12, multirev_long 9.5e-13 / 1.3e-12, hyperbolic_fast 2.8e-15 / 3.0e-15. All checks pass: standalone 3/3 (OpenMP on and off), Lamberts-Comparison 9/9, MATLAB suites (hard mixed fixture 3,004/3,009), identical statuses to commit 36228ca on 400,000 extreme requests and on the 1M hard mixed dataset.

**Floors set by the double-precision inputs.** Solving a case exactly (quad precision) from its rounded inputs gives the smallest error any double-precision solver can guarantee against the generating orbit. For the multirev_long case that sets EXCALI's maximum (m = 4, 7.7e-11 above the minimum time, which ivLam fails), that floor is 4.0e-12; EXCALI gives 1.75e-11 and confluentcontour 8e-13 (below the floor, by cancellation of rounding errors). ivLam's 1.28e-12 maximum in that class is over the cases it solves, so no solver can beat it there except by chance. For the stress cases near_collinear and near_half_long the floor is 1.6e-9 and 3e-10, above the 1e-10 velocity tolerance; singular_half and singular_coincident have no defined transfer plane. Medians are at the rounding floor for all solvers (exact-solution floors 1.2-1.4e-16; ivLam 2.0-2.7e-16; EXCALI 2.0-4.5e-16); EXCALI's double-double f,g velocities reach about 2.0-2.4e-16 in every class but cost about 13 ns per case if used everywhere, so they are used only near the parabola and the short-way branch point.

## Smaller starter tables — 2026-10-27 (later)

The zero-height rectangle needs a starter only within a few percent of the singularity distance, so the starter grids were made coarser: 32 points per axis for the zero-revolution tables (64 before) and 16 × 16 for the multi-revolution branch roots (32 × 32). `src/excali_tables.inc` shrinks from 9,035 to 2,735 lines. The solver code is unchanged.

The multi-revolution minimum-time table is not a starter: it decides feasibility outside a 1e-6 margin. Measured against bisection truth for m = 1–40 and sampled m up to 1000, the previous 256-point table's worst relative error was 4.4e-7, at m = 5 (interpolated between the m = 4 and m = 6 slices), not "below 1e-7" as recorded earlier; it held below 1e-7 only on the slices. That was still inside the margin. The table now has an m = 5 slice and 128 points, with a worst error of 1.04e-7. With 64 points it exceeds the margin (1.2e-6), and with 16 it gave 108 wrong infeasible statuses on the mixed set.

Coarser options measured (single-thread A/B against the previous tables on the same cases, Apple M1 Pro; an identical build measures 0.999–1.022, so differences under about 2% are noise):

| Tables (zero-rev / branch roots / t_min) | Lines | Mixed set, single / batch | Hard 1M set, single | Statuses |
|---|---:|---|---|---|
| 32 / 16 / 128 + m = 5 slice (shipped) | 2,735 | 1.00 / 1.00 | 1.03 | identical |
| 16 / 16 / 128 | 1,641 | 1.00–1.03 / 1.04 | 1.08–1.12 (hyperbolas 1.19–1.26) | identical |
| 8 / 8 / 16 | 359 | 1.19 / 1.31 | — | 108 wrong |

8-node rectangles of half-width 1/8 of the bound everywhere, so that a coarser starter still finishes in one rectangle, were not faster. With the shipped tables the one-rectangle path finishes 99.4–100% of each mixed class (20,000 cases each), as before; about 1% of zero-revolution elliptic cases with `v >= 0.3` now miss the 3-node window and take the general path.

**Checks:** identical statuses to the previous tables on the 220,000-case mixed set and the 1,000,009-case hard mixed dataset (scalar and batch), with solutions within 8e-13 of each other; the same 12 hard-set cases fail with both.

## No Halley search: fixed-node multi-revolution formulas (2026-11-03)

The Halley search for the multi-revolution minimum time was removed. Where the tables do not decide (within the 1e-6 margin of `t_min`, outside their `tau` range, or a rejected rectangle), two fixed 32-node zero-height ratios run instead: one for the minimum location `v*`, as the zero of `T'(v) sin^4 v`, and one for the branch, in `rho = sign(v-v*) sqrt(1-(t_min/T)^(1/3))`. Both use `[0, pi]`, or a scale-aware interval around an asymptotic estimate when `|b/(y0+b)| > 0.95`. The method and the experiments behind it are in Lamberts-Comparison `research/closed-form/`.

- **Path counts.** On the 119,531 truth cases of Lamberts-Comparison `research/negrete-reports/data` (zero_rev, adversarial, multirev, infeasible, harness_mixed and dev) (all multi-revolution, near-`t_min`, high-`n`, infeasible and adversarial sets), the Halley code is never reached once the asymptotic interval is used. On a 24,000-case stress set (transfer angles down to 1e-12 from 0 and 2 pi, radius ratios to 0.01, m = 1 to 1000, times from 1e-12 above `t_min` to 100 times it, plus infeasible times), the certificate never failed. The 780 adversarial near-`t_min` cases are now solved with no bracketed contour; it took 400 of them before.
- **Accuracy.** Near-`t_min` adversarial cases: 172 above 1e-10 (195 before), worst 1.1e-8 (2.6e-8 before); these are at the conditioning floor. All other categories are unchanged. On the stress set, the previous version failed 12 cases (status -4) and the new one none. Where velocities differed by more than 1e-8 (578 cases), quad-precision Kepler propagation put the new answer 10x closer to `r2` in 155 cases and the old one in 17; the worst arrival errors are 2.9e-9 and 8.4e-9.
- **Tests.** Lamberts-Comparison CTest 9/9 (Release and Debug, including the 32,400 feasibility decisions), standalone 3/3.
- **benchmark_mixed** (2,000 cases per class, M1 Pro, two interleaved runs): zero failures in every class. Overall 0.87-0.88 of ivLam before and after. Multi-revolution classes take 1-2% longer, within run-to-run noise. Maximum errors are equal in 8 classes, lower in `multirev_1` (5.1e-13 vs 5.8e-13) and `multirev_long` (8.2e-13 vs 9.0e-13), and higher in `multirev_2_3` (2.0e-13 vs 1.0e-13).
- **Cost.** The adversarial near-`t_min` class takes 864 ns per solve against 365 ns with Halley: 32 time evaluations for `v*` cost more than one or two Halley steps. Without the multi-revolution tables, every multi-revolution solve would take this path (1,159 ns against 204 ns on random cases), so the tables stay.
- **Still not iteration-free.** Zero-revolution starter misses at `b/(y0+b) -> 1` (near-parabolic orbits, very small transfer angles) and very short hyperbolas still use the bracketed contour. One stress case of multi-revolution branches far from `t_min` at `|b/(y0+b)| > 0.95` reached it too.
- **Size.** `src/excali_solver.f90` grew from 1,773 to 1,939 lines.

## General path replaced by fixed formulas (2026-11-04)

The general path is replaced by fixed formulas. These are deleted:
- the closed-form zero-revolution starters;
- the 5-pass stage control (`rect_root`, `canonical`, `reach`, `contour_root`, `pole_coordinate`);
- both bracketed contours (`converge` in the anomaly, `converge_y` in `sqrt(y)`) with `ytime` and `xof`;
- the multi-revolution branch starter (`branch_start`) and the helpers only they used (`hermite`, `rasinh`, `rasin`, `qlog`).

Every case the tables do not settle now goes to one fixed 32-point zero-height ratio chosen by region (README, General path). `src/excali_solver.f90` went from 1,773 lines (before 2026-11-03) and 1,939 (after the Halley removal) to 1,573.

- **Fallback population.** The 720 zero-revolution cases of the truth sets that previously reached a bracketed contour are all solved by the fixed formulas:
  - 715 short-way: near-parabolic orbits, `b/(y0+b) -> 1`, near-collinear, very short hyperbolas, in `Y = sqrt(y)`;
  - 5 long-way at `b/(y0+b) ~ -0.99997`, in `zeta`.
- **Truth sets** (119,531 cases, Lamberts-Comparison `research/negrete-reports/data`): no hard failures. Every category has the same or fewer failures above 1e-10 than the original:

  | Category | Original (fails above 1e-10 / max error) | Fixed formulas |
  |---|---|---|
  | near-parabolic orbits | 6 / 1.0e-8 | 5 / 6.1e-10 |
  | adversarial `lambda -> +1` | 3 / 1.8e-9 | 0 / 5.8e-11 |
  | adversarial near-collinear | 243 / 5.4e-4 | 240 / 1.3e-6 |
  | adversarial near-`t_min` | 195 / 2.6e-8 | 172 / 1.1e-8 |

  The other categories are equal within rounding; the largest change is 6.1e-14 against 2.8e-14 in one dev class.
- **Zero-revolution stress set** (40,000 cases: transfer angles within 1e-12 of 0, pi and 2 pi, radius ratios to 0.01, both directions, `t/t_p` from 1e-6 to 1e9):
  - 0 failures (the original: 1).
  - Where velocities differ by more than 1e-10 (527 cases), quad-precision Kepler propagation puts the new answer 10x closer to `r2` in all 527.
  - At `t/t_p > 1e5` both versions have the same, input-limited, arrival error.
- **Multi-revolution stress set** (24,000 cases): 0 failures (the original: 12).
- **Tests.** Lamberts-Comparison CTest 9/9 (Release and Debug), standalone 3/3.
- **benchmark_mixed** (2,000 per class, two interleaved runs): zero failures. Overall 0.86-0.88 of ivLam against 0.86 before, within run-to-run noise. Maximum errors are equal in 8 classes and lower in `multirev_1` and `multirev_long`; `multirev_2_3` is 2.0e-13 against 1.0e-13.
- **Cost.** A case on the general path now evaluates `T` at 32 nodes, plus at most 6 more for the region choice and the `zeta` piece selection.
  - Random truth classes (where the tables decide almost everything): 3-5% slower.
  - Hard categories, which mostly take the general path:

    | Category | Before | After |
    |---|---:|---:|
    | near-parabolic orbits | 262 ns | 540 ns |
    | `lambda -> +1` | 361 ns | 755 ns |
    | `lambda -> -1` | 211 ns | 1,192 ns |
    | near-collinear | 578 ns | 883 ns |
    | near-`t_min` | 390 ns | 868 ns |

    These are single runs, M1 Pro under load.

## Single path for zero-revolution ellipses (2026-11-04)

Zero-revolution ellipses take one path: geometry → starter → one fixed Gauss–Legendre rectangle → `Q1/Q0` → velocities. The rectangle's result is final: it has no certificate, retry or fallback. Ellipses outside the tables (`|tau| > 0.9`) go directly to their region formula.

- **Normalized elliptic table.** The table now stores `(pi-v)/(s*ck)` and `v^2*T0'/(t-t_p)`, the root divided by its pole and parabola asymptotes. Its error relative to the singularity bound `min(v, pi-v)` falls from 10% to 3.9% (worst case up to `|tau| = 0.9`), and to 0.11% for `|tau| < 0.75`. The table size is unchanged (32×32×2), and the hyperbolic and multi-revolution tables are byte-identical.
- **Node count, fixed before evaluation.** 3 nodes at half-width bound/256 for `|tau| < 0.75`; 8 nodes at bound/8 for `|tau| >= 0.75` and the near-parabolic series. The near-parabolic window widens from 3% to 3.1% of `t_p`, so it meets the table at `s = 0.99`.
- **Measured starter errors** (66,114 ellipses of the truth sets, relative to `min(v, pi-v)`):

  | Starter | Median | p99 | Max |
  |---|---:|---:|---:|
  | Previous table | 3e-5 | 4.2e-2 | 1.0e-1 |
  | Near-parabolic series | — | — | 7.8e-2 |
  | Closed-form Hermite | — | — | 6.1 |

  The closed-form Hermite starter of earlier versions fails near `|b/(y0+b)| -> 1`.
- **Alternatives tried:**
  - Without the normalized table, one 8-node rectangle at bound/8 misses roots (34 failures in 40,000). At bound/4 it solves everything but loses accuracy (worst 3.5e-12 against 1.2e-13).
  - A 3-node rule everywhere needs a starter within 0.16%.
  - Tabulating the two scale factors in `tau` was too inaccurate: `t_p/T0'` spans five orders of magnitude and interpolates to 21% error.
  - The kept form computes the one factor its side of the table needs.
- **Truth sets** (119,531 cases): failure counts and maximum errors identical to the previous version in every category.
- **Stress sets:** 40,000 zero-revolution and 24,000 multi-revolution cases, 0 failures (the original: 1 and 12).
- **Tests.** Lamberts-Comparison CTest 9/9 (Release, Debug), standalone 3/3.
- **benchmark_mixed** (3 interleaved runs, against the original EXCALI before 2026-11-03): zero failures. Elliptic classes take 0.99–1.035 of the original's time, within this machine's run-to-run spread. A controlled circular-transfer timing gives 153 ns for both.
- **Code generation.** Adding the 8-node branch inside `fast_coordinate` cost 14 ns on circular transfers that never take it, so `solve_chunk` widens the rectangle after the call instead.
- **24-point region formulas.** The per-region formulas (`ratio_root`) use 24 Gauss–Legendre points instead of 32.
  - Accuracy is identical in every truth category; 16 points lose accuracy (up to 4.9e-12).
  - The multi-revolution minimum keeps 32 points.
  - `excali_batch` on the 70,000 elliptic truth cases (zero_rev categories 1, 3, 4; one thread):

    | Version | Time per case |
    |---|---:|
    | Original EXCALI | 140 ns |
    | Previous version | 207 ns |
    | Single path, 32 points | 200 ns |
    | Single path, 24 points | 186 ns |

    The gap to the original comes from the degenerate ellipses of category 3 (chord under 1% of the semiperimeter, `|tau| > 0.9`). They use the region formulas, where the original used starters, rectangles and the bracketed contour.
- **Single-solve profile** of a circular transfer: geometry and table coordinates 37 ns; table lookup, coordinate and centre trig 50 ns; rectangle 42 ns; acceptance 15 ns; velocities 10 ns.
- **Measured and not kept** (circular-transfer timing, 153 ns):

  | Change | Saving |
  |---|---|
  | Skipping the second moment and error estimate for ellipses | 1 ns or less |
  | A three-term shift series in `rotate` | 0.5 ns |
  | A specialized, unrolled 3-node elliptic rectangle | 2 ns, for about 30 more lines |
  | Tabulated table scales | made the table up to 21% wrong |
- **Short-way table toward `tau -> 1` (`tabp`).** A second 32×32×2 table covers short-way ellipses for `tau` in `[0.75, 1-1e-6]`, on a grid uniform in `1-sqrt(1-tau)`. It splits at `v = pi/2` (`T_half = T(pi/2)`, closed form). On the pole side it stores `(pi-v)/(s1*(K/T_half)^(1/3))`; `(K/T_half)^(1/3)` stays finite as `tau -> 1`. On the parabola side it stores `v^2/(s2*(s2+2*e0))`, `s2 = (t-t_p)/(T_half-t_p)`, `e0 = t_p/(T_half-t_p)`. As `lam -> 1`, `v^2` is proportional to `s2*(s2+2*e0)`, which absorbs the branch point's boundary layer.
  - Starter error: at most 3.1e-5 of `min(v, pi-v)` over the whole range, so these ellipses use 3 nodes.
  - Truth sets: identical or better in every category. Stress sets: 0 failures.
  - `excali_batch` on the 70,000 hard elliptic truth cases: 187 → 147 ns per case (the original EXCALI: 140 ns). The short-way degenerate ellipses that used the 24-point Y formula now take the single path.
- **Long-way closed-form starter (`tau < -0.9`).** As `tau -> -1`, `T/|b|^1.5` has a plateau at `pi/4` between a body, where `d = pi - v` depends only on `th = t/|b|^1.5`, and a pole layer of width `eta = sqrt(2*y_end/|b|)`. `long_d` uses three forms:
  - the pole layer `(pi/4)((1+eta^2/d^2)^1.5-1) = th - pi/4`, solved exactly when it gives `d < 2*eta`;
  - the transition cubic `d^3/2 + (th-pi/4) d^2 = 3*pi*eta^2/8`, in the scaled form `x^3/2 + chi x^2 = 3*pi/8` (`d = eta^(2/3) x`), solved through `1/x` for `chi > 0` and directly otherwise by trigonometric or hyperbolic formulas (cancellation-free), when it gives `d < 0.1`;
  - otherwise the body limit from a 64-point table (`tabb`).

  Results:
  - Starter error at most 6.4% of `min(v, pi-v)` for `tau` from -0.9 to -0.999999; these cases use 8 nodes at bound/4.
  - Short-way `tau > 1-1e-6` clamps `tabp` at its edge, also 8 nodes at bound/4.
  - No zero-revolution ellipse in the truth sets or the 40,000-case stress set reaches the region formulas, so their elliptic dispatch was removed.
  - Truth sets equal or better in every category, both stress sets 0 failures, all tests pass.
  - `excali_batch` on the hard elliptic set: 145 ns per case (original EXCALI 140 ns).
  - Earlier attempts that failed: the `v = pi/2` split on the long way (errors to 1e6), and a matched composite with a correction table (Cardano cancellation deep in the layer).
- **benchmark_mixed** (final, 3 interleaved runs): zero failures; maximum errors equal or lower in every elliptic class. Overall time is 1.02 of the original (174 against 171 ns); elliptic classes 1.02–1.07. The controlled circular-transfer timing is equal (153 ns), so part of this is code layout between separately built binaries.
- **Size.** `src/excali_solver.f90` is 1,687 lines, against 1,573 before. The additions are the two table lookups (`table_vp`, `long_d` with `depressed_max`), the 24-point constants, the per-side table scale and the node widening. The removals are the elliptic region dispatch, now unreachable.

## Single path for zero-revolution hyperbolas (2026-11-04)

Every zero-revolution case, elliptic or hyperbolic, now takes the same path: starter, one rectangle, `Q1/Q0`, velocities, with no certificate. The zero-revolution region formulas (`fixed_root`, the `Y`, `X` and `w` kinds) are removed; `general_root` serves multiple revolutions only.

- **Where hyperbolas left the table path before** (logging build, `e290b01`):

  | Set | Hyperbolas to region formulas | Cause |
  |---|---:|---|
  | zero_rev (100,000) | 2,091 | 2,010 short-way `tau > 0.9` (off `tabh`); 25 `abs(tau) < 4e-4`; 56 certificate rejections near the parabola (error estimate 4–60 ulp of a tiny `X`) |
  | adversarial (4,771) | 887 | off the table at `abs(tau) > 0.9` (559) or `abs(tau) < 4e-4` (151), certificate rejections (177) |

  The old `xi` table (`u < 0.1`) was also wrong by up to 4.7 times the singularity distance for short-way transfers within about 0.3 degree of 180 degrees; the certificate caught those.
- **New hyperbolic table** (`tabh`, 32×32×2, replaces the 32×32×3 one): `u = t/t_p` in `[0.1, 1]`, `sqrt|tau|` in `[0, sqrt(1-1e-6)]`, clamped beyond. Short way `w^2/((1-u) wn^2)` with `wn^2 = L^2 (1 - exp(-wb^2/L^2))`, `L = 1 - log u`; long way `w^2/((1-u) L^2)`. Starter error at most 6.3e-4 of the singularity bound for `tau` from 1e-16 to 1 - 1e-9 (the rational blend `wb^2 L^2/(wb^2 + L^2)` reaches 1.3% near `tau = 4e-5`, because `wb ~ log(1/tau)`). These hyperbolas use 3 nodes at bound/256.
- **Closed form for `u < 0.1`.** `T = sqrt2 c sqrt(y)/E`, `E = e^w`, `c = y0 + b`, `y = c - bE/2` (the large-`w` estimate the region formula used to centre its interval). Measured starter error, relative to the singularity bound:

  | Branch | Coordinate | Worst error |
  |---|---|---:|
  | long way | `w = log E` | 1.1e-2 (at `u = 0.1`; 2.3e-4 at `u = 0.01`) |
  | short way, `y >= y0/2` | `w = log E` | 1.1e-2 |
  | short way, `y < y0/2` | `xi` from `y = (tE/(sqrt2 c))^2` via `cosh w = (c-y)/b` | 1.1e-2 |

  These use 8 nodes at bound/8.
- **Bugs found on the way:**
  - The short-way bound `min(2w, wb-w)` had no cap. At large `w`, `T ~ exp(-w)` has near-copies of the root at `w ± 2 pi i`, and the half-width could exceed 1, past the node tables (`kr < 0`). The bound is now `min(2w, wb-w, pi)`, matching the long way's `min(w, pi)`.
  - `xi ~ t^2` reaches half-widths below `2^-62` (the node tables' range) at `t/t_p ~ 1e-9`; `rect_ratio` now uses `cos = cosh = 1`, `sin = sinh = L u/2` there, exact to rounding.
- **Rejected:** the `Y = sqrt(y)` coordinate for the short-way closed form. It is accurate (no other zero within 0.98 of `sqrt(y(pi)) - Y`, checked by the argument principle for `tau` from 1e-8 to 1), but needs an `asinh` per node: the fast-hyperbola category took 598 ns against 234. `xi` for every short-way closed form fails where `xi` is large: `sinh w` by angle addition cancels from about `1e20` (errors up to 179 in the stress set).
- **Truth sets** (119,531 cases): failure counts equal or lower in every category, maximum errors equal within rounding (largest change 1.74e-12 → 1.90e-12 in zero_rev category 2); per case, none more than 4 times worse above 1e-13, two better.
- **Stress sets** (quad-precision Kepler judge): the 40,000-case zero-revolution set and a new 100,000-case hyperbolic set (`t/t_p` from 1e-9 to 1, `r2/r1` from 1e-3 to 1, transfer angles clustered within 1e-12 of 0, 180 and 360 degrees): 0 failures, no velocity differing from the previous version by more than 1e-10.
- **Tests.** Lamberts-Comparison CTest 9/9 (Release, Debug), standalone 3/3.
- **Speed** (one thread):

  | Measure | Previous | Now |
  |---|---:|---:|
  | `excali_batch`, adversarial set | 404–414 ns | 259–262 ns |
  | `excali_batch`, zero_rev set | 152–155 ns | 140–143 ns |
  | `excali_batch`, multirev / harness_mixed | 170 / 123–125 ns | 170 / 124–128 ns |
  | single solves, adversarial hyperbolic categories 11–14, 18 | 356–483 ns | 182–225 ns |
  | single solves, fastest hyperbolas (category 15) | 234 ns | 265 ns |
  | benchmark_mixed, every class | — | within run-to-run spread of the previous version |

- **Size.** `src/excali_solver.f90` is 1,608 lines (1,687 before): `fixed_root`, `xlow`, three `zeval` kinds, `table_xi`, `tabh_at` and the `yy` velocity path are gone; the closed form adds 22 lines.

## One region selector for zero revolutions (2026-11-04)

`zero_rev_start` replaces `table_coordinates`, `table_root` and the node-widening rules in `solve_chunk`. It returns the coordinate, the starter and an accuracy class (0: 3 nodes at bound/256; 1: 8 nodes at bound/8; 2: 8 nodes at bound/4), and `fast_coordinate` alone turns the class into the rectangle. Every zero-revolution region, its starter and its rectangle are listed in one comment table (and in the README).

- **Results identical** in 119,526 of 119,531 truth cases (bit for bit). The other 5 (adversarial categories 11 and 14) are near-parabolic hyperbolas with `w` just above 0.3. They now get the 8-node rectangle the series starter is meant to have; before, they fell through to 3 nodes at bound/256. Both versions are at rounding (2e-15 to 2e-14).
- **Stress sets** (40,000 zero-revolution, 100,000 hyperbolic): 0 failures, no velocity differing from the previous version by more than 1e-10. Tests: CTest 9/9 (Release, Debug), standalone 3/3.
- **Code generation.** Calling the selector inside the geometry loop cost 8–13 ns per solve. A separate loop costs nothing measurable. Splitting the hyperbolic closed form into its own routine, a module-level selector, or arithmetic widths instead of `select case` changed nothing.
- **Speed** (best of 3–15 interleaved runs, against the previous commit):
  - `excali_batch`: equal or faster on every set (zero_rev 139 → 137 ns, adversarial 256 → 248, multirev 169 → 168).
  - Single solves: equal or 1–8% faster in every truth category.
  - Circular-transfer microbenchmark: best 152 → 157 ns, mean 158 → 162, within this machine's run-to-run spread (the same binary measured 152–166).
  - benchmark_mixed: every class within 1 ns of the previous version, same maximum errors.
- **Size.** `src/excali_solver.f90`: 1,601 lines (1,608).

## Multiple revolutions on a single path inside the tables (2026-11-04)

The multi-revolution table rectangle no longer carries a certificate, and the certificate code is gone (`rect_ratio` computes `Q0` and `Q1` only; `fast_accept` only moves the root). `general_root` now serves multiple revolutions off the tables only (`abs(tau) > 0.9`): 28 cases of the truth sets, against 639 before.

- **Normalized table.** The old `mr_v` (16×16 per slice, `v` itself) was off by up to 8 times `abs(v - v*)` near the fold and 20% of `min(v*, pi - v*)` elsewhere; the certificate hid that. The new `mr_q(32,16,2,11)` stores `q = (v - v*)/(sigma (P - v*))`, `P` the pole on the branch's side, on one half-grid in `abs(sigma)` per branch.
  - **Normalizations tried** (worst starter error relative to the rectangle's bound, 1,500 random cases including `m` between slices up to 1,000):

    | Table | Near the fold (`1e-3 <= abs(sigma) < 1e-2`) | `abs(sigma) > 0.9` |
    |---|---:|---:|
    | `v`, 16×16 (old) | 800% | 0.5% |
    | `(v - v*)/sigma`, 32×32 | 2.2% | 15% (interpolated `v*` error at the pole) |
    | `q`, one grid across `sigma = 0` | 1,500% (`q` jumps at `sigma = 0`) | 0.15% |
    | `q`, half-grids, 16 `tau` | 7% | 0.3% |
    | `q`, half-grids, 32 `tau` (kept) | 2.5% | 0.17% |

  - End-to-end in Python (3,000 cases, `abs(sigma) >= 1e-2`): one 8-node rectangle at bound/8 gives every root to the truth's precision (3e-13); the starter is at most 0.49 of the half-width from the root.
- **Near the fold** (`abs(sigma) < 1e-2`, or `t` within the 1e-6 margin of the table's `t_min`): `mrev_fold`.
  - One 8-node ratio of `D = T' sin^4 v` at `min(v*, pi - v*)/8` about the table's `v*` (within 6e-5 of that bound) gives `v*`, and `t_min = T(v*)` decides feasibility exactly.
  - The starter is the fold series `v* + k s (1 - T''' k s/(6 T''))`; in Python it is within 0.6% of `abs(v - v*)` for `abs(s) < 3e-2`.
  - It replaces the 32 + 24-node `mrev_global`/`branch_fixed` pair for these cases.
- **Rejected:**
  - A table-free starter, a cubic Hermite inverse of `rho(v)` from the pole and fold slopes. It fails for `abs(tau) > 0.7`, where the `sqrt(y)` branch points crowd an end.
  - The two-root rectangle (`research/closed-form/two-root` in Lamberts-Comparison), which only contains both roots for `abs(sigma) < 1e-2`.
- **Truth sets:**
  - multirev: identical summary (0 failures, maximum 3.66e-12); per case 9 better and 3 worse by more than 4 times (worst 8e-13).
  - Near-minimum adversarial category 10: failures above 1e-10 drop from 172 to 161; the maximum goes from 1.1e-8 to 1.6e-8; per case 83 better, 76 worse. These are the cases within a few ulp of `t_min`, where `abs(v - v*)` is about `sqrt(eps)` and both methods are limited by the rounding of `t_min` itself.
  - Every zero-revolution result is bit-identical to the previous commit except one (one rounding, from computing `Q1/Q0` directly).
- **Stress sets:**
  - Multi-revolution, 24,000 cases including 6,000 infeasible: 0 failures, 0 feasibility disagreements. Where velocities differ by more than 1e-8 (483 cases, all near the fold), the new version is 10 times closer to `r2` in 44 and the previous one in 29; worst arrival error 2.6e-9 against 2.9e-9.
  - Zero-revolution (40,000) and hyperbolic (100,000): unchanged, 0 failures.
- **Tests:** CTest 9/9 (Release, Debug), standalone 3/3. benchmark_mixed maximum errors: multirev_1 5.1e-13 → 7.9e-14, multirev_long 8.2e-13 → 2.8e-13, multirev_2_3 2.0e-13 → 2.3e-13.
- **Speed** (best of 3, one thread):

  | Measure | Before | Now |
  |---|---:|---:|
  | Single solves, near-minimum category 10 | 899 ns | 456 ns |
  | Single solves, multirev set | 214 ns | 215 ns |
  | `excali_batch`, multirev set | 169 ns | 163 ns |
  | `excali_batch`, adversarial set | 257 ns | 165 ns |
  | Circular-transfer microbenchmark, mean of 15 (best) | 162 ns (157) | 163 ns (162) |

  Making `mrev_fold` a module procedure instead of an internal one was slower (166 ns).
- **Size:** `src/excali_tables.inc` grows by about 1,860 lines (the 11,264-value `mr_q` against 2,816 for `mr_v`).

## Near-parabolic series final at the rounding floor (2026-11-05)

Lamberts-Comparison CI (Linux x86-64) failed `test_excali`: 4 of 10,011 cases within 1e-12 of the parabolic time returned status -4. With `X` around 1e-13, `T - t` at the rectangle's nodes is a few hundred ulps, and on x86-64 one node rounded to exactly 0, so `Q1/Q0` was NaN. The region formula used to catch these; since the single path there was no recovery.

- **Fix.** `near_parabola_x` also returns the relative size `rc` of its second-order term. When `|rc| < 1e-8` the dropped third-order term is about `rc^2`, so the series is the root to rounding and no rectangle is evaluated (accuracy class -1).
- **First attempt rejected:** an absolute threshold `|X| < 1e-6` raised zero_rev category 3 (degenerate ellipses) from 5 to 115 cases above 1e-10, because there `T''/T'` is large and the series is not converged at that `X`.
- **Near-parabolic stress set** (60,000 cases, `t/t_p - 1 = ±1e-16..1e-2`, every geometry of the zero-revolution stress set):
  - previous version: 0 failures on macOS, 1,551 on Linux x86-64 (gcc 14, Docker amd64);
  - this version: 0 on both.
  - Maximum arrival error 1.8e-12 on macOS (3.5e-12 before) and 1.7e-12 on x86-64 (4.3e-12).
- **test_excali on x86-64:** parabolic worst 3.6e-14, 0 bad. Truth sets unchanged in every category.

## Multi-revolution robustness and the dead whole-interval path (2026-11-05)

A larger multi-revolution stress run (120,000 cases), plus a run with times up to 1e48 `t_min`, found three faults:

1. **Region-formula failures.** 5 cases failed (status -4) in both this version and the previous one. All were right-branch roots at transfer angles within 1e-3 of 0 or 360 degrees. `branch_fixed`'s one 24-node ratio over the whole `d` interval (length about 3) ended about 1e-9 short of its own test, with the root near the end close to the `v = 0` pole. The `d` and `v` intervals are now cut into pieces of length at most 1/2, as `zeta` already was; error estimates fall to 3e-15 or less and the failures to 0.
2. **Huge times.** At `t` up to 1e48 `t_min`, 3,116 of 18,000 cases fell back to `general_root`. The starter `v* + |sigma| (P - v*) q` rounded onto the pole. The table now stores `g = (q - 1)/(1 - |sigma|)`, which is finite at the pole, and the starter is formed as the pole distance `(1 - s)(1 - s g)|P - v*|` with `1 - s = r/(1 + s)`. That is exact as `t -> infinity`, and the starter comes back in the rectangle's own coordinate. Storing `q - 1` alone still left a 1e-9 offset at the pole (from the grid's `1 - 1e-9` end sample).
3. **The fold boundary.** At `t = t_min (1 - 16 eps)` the infeasibility test and the fold test rounded in opposite directions, and `sqrt(1 - t_min/t)` was NaN. Both tests now use the single quantity `t_min/t`.

**Dead code.** `general_root` is reached only off the multi-revolution tables (`abs(tau) > 0.9`). There `abs(b/(y0 + b)) >= 0.999985` over 4 million sampled geometries, so the whole-interval branch of `mrev_global` (`abs(b/(y0+b)) <= 0.95`) and its node tables are gone.

**Results:**
- Huge times: 0 fallbacks from `t_min` to 1e48 `t_min`.
- Multi-revolution stress: 0 failures in 120,000 cases (the previous version: 5) and 0 in 24,000.
- Zero-revolution, hyperbolic and near-parabolic stress sets: unchanged, 0 failures.
- Truth sets: unchanged in every category.
- `test_excali` on Linux x86-64 (Docker): 0 bad.
- CTest: 9/9 and 3/3 (Release, Debug).
- Batch timings: equal to the previous version.

## Confidence flag on the final rectangle (2026-11-05)

Every final rectangle again computes its second moment. `de = |Q2/Q0 - (Q1/Q0)^2|` vanishes for a single pole. When `de/L > 1e-8 R` (`L` the half-width, `R` the singularity bound), the case returns status 2, "solved, low confidence". The velocities are returned and nothing is re-solved. It costs about 1.5 ns per single solve (circular microbenchmark, best of 15) and 3-5 ns per case in batch.

- **Where it appears.** It never appears in the zero-revolution truth sets (largest `de/(R L)` = 4.8e-11), harness_mixed, dev, multirev (1.4e-9), the 40,000 zero-revolution, 100,000 hyperbolic or 60,000 near-parabolic stress sets, or in CTest.
  - It appears only within about 3e-9 of a multi-revolution minimum time: 309 cases of adversarial category 10, which places `t` a few ulps from `t_min`, and stress cases with `t/t_min - 1` from 1e-12 to 2.6e-9.
  - There the root is `sqrt(t - t_min)`-conditioned. The flagged cases' worst error is 1.6e-8, against 2.7e-11 for the unflagged ones in the same category, so the flag marks exactly the inherently inaccurate results.
- **Tests.** Status 2 counts as a failure in every test outside that zone:
  - `test_excali_standalone` and the Lamberts-Comparison tests accept only status 0 or 1;
  - `test_excali_minimum` checks `t/t_min - 1` from 1e-6 to 1e-1;
  - the stress harnesses count any other status as a failure.

  A future starter or width regression therefore shows up as a test failure, not as a wrong velocity. MATLAB's accuracy tables count status 2 as a valid solution and measure its error as usual.

## Single path for multiple revolutions off the tables (2026-11-05)

`mrev_edge` replaces `general_root` (`mrev_global`, `branch_fixed`, `ratio_root`, `zeval`) and the 24- and 32-point rules. Multiple revolutions off the tables (`abs(tau) > 0.9`) now take the same three steps as everything else, and every case in EXCALI takes one path.

- **Coordinate.** `zeta` about the crowded pole, `sin(x/2) = s0 sinh(zeta)`, `x = v` (short way) or `pi - v` (long way). There `y = y_c cosh^2(zeta)` exactly, and `T = A coth^3(zeta)/cos^3(x/2) + sqrt(y_c) cosh(zeta) (±y_c cosh^2(zeta) G(x) + a)`. Each term is evaluated directly, with `zeta`-derivatives to third order by Taylor arithmetic (`jmul`, `jcomp`). In `v`, `T'` near `v*` cancels about 30 digits once `1 - abs(tau)` falls below about `1e-10`, so a double-precision formula in `v` cannot locate the minimum there.
- **Minimum time.** One 8-node ratio of `1/T'` over half-width 2 about the large-`m` estimate puts `v*` within 3.7e-5 of its bound in `zeta` (1.8e-9 in `x`). A second 8-node ratio at bound/16 then reproduces a 32-node ratio's `v*` to rounding (3.7e-16 relative). Feasibility decided by `t_min = T(v*)` matched 50-digit minima within 2e-15 to 4e-15 in 40 random cases. That is inside the 16-ulp `t = t_min` tolerance.
- **Starter** (`Delta = (t - t_min)/abs(t_min - A)`, `A = n pi abs(b)^1.5/4`):
  - Near the fold: the fold series in `zeta`.
  - Crowded side, short way: the layer cubic `2w^3 - D w^2 + 1 = 0` fitted at `zeta*`. Below `zeta = 2.5` it uses `alpha coth^3(zeta) + beta cosh(zeta)` fitted exactly at `zeta*`, with `cosh` taken at the layer-only root. The plain `atanh((A/t)^(1/3))` reached 0.24 of the bound at `abs(tau) = 0.91`; the fit cuts that to 0.04.
  - Crowded side, long way: `atanh`, the fold model, the plateau cubic or the anchored body.
  - Far side: the body limits simplify to `4 theta cos^3(x/2) = n pi + (x + sin x)` (short way) and `n pi - (x + sin x)` (long way).
    - Short way: three fixed Newton steps in log form invert it (6e-6; replaced by a table on 2026-11-06, see below). Within `8 v*` of the far pole the inverse is final; otherwise its secant slope from `v*` sets the fold cubic. The cubic alone overshoots by about `v*/2` as the root nears the pole: 0.34 to 3.5 of the bound for `(pi - v)/v* < 2`, where the body inverse is within 1e-3.
    - Long way: the body minimum takes two Newton steps (5e-8). Each branch takes two steps in `sqrt(log Tb - log Tb_min)`, which is regular through that minimum (1e-6). (Replaced 2026-11-06, see below.)
  - The four-term fold model takes three steps in `x` or `1/x` (1e-4). (Replaced by Ferrari's formula 2026-11-06.)
  - All cubics are solved in closed form; `cubic_pos` also replaces `long_d`'s inline cubic.
- **Rectangle.** 8 nodes at bound/4. It is in `zeta` when the starter has `x < 1`, otherwise in `x`, or in `pi - x` above `3 pi/4`. In `x` near `pi`, `pi - x` keeps only absolute precision: the far-branch velocities were up to 30 times less accurate than the old path's until this was added.

Validation:

- **New stress set, `abs(tau) > 0.9`.** 10,000 geometries × 4 time bands:
  - transfer angle `10^-7` to `0.4` rad from 0 or 360 degrees, `1 - r2/r1` from `1e-9` to `0.1` or 0;
  - both directions (an earlier run used direction +1 only, which never reaches the long way);
  - `m` from 1 to 30,000;
  - `t/t_min - 1` in `[1e-12, 1e-8]`, `[1e-8, 1e-4]` and `[1e-4, 1e4]`, plus infeasible times.

  | | old path | `mrev_edge` |
  |---|---|---|
  | status -4 | 5,890 of 40,000 | 0 |
  | status 2 | 0 | 6,936, all with `t/t_min - 1 <= 2.7e-9` |
  | worst arrival error, `t/t_min - 1` in `[1e-12, 1e-8]` | 5.4e-12 | 5.6e-12 |
  | worst arrival error, `[1e-8, 1e-4]` | 5.2e-9 | 6.1e-12 |
  | worst arrival error, `[1e-4, 1e4]` | 3.4e-4 | 3.4e-4 (same case) |
  | 10 times worse than the other, where both solve | 205 | 38 |

  Arrival errors are by quad-precision Kepler propagation. In all 38 cases where `mrev_edge` was 10 times worse, the two velocity vectors differ by at most 4.7e-15. Arrival errors at large `m` and long `t` sit at the rounding floor: in the worst case, one ulp in `v1` moves the arrival error from 1.0e-4 to 2.75e-4. Feasibility disagreed in 693 cases. All were below `t_min`, where the old path failed with -4 and `mrev_edge` returned 5. In 29 sampled cases, 50-digit minima put `t` below `t_min` by 5e-9 to 1e-5.
- **Extreme set.** 5,000 geometries with transfer angles down to `1e-12` rad and `1 - r2/r1` down to `1e-15`. The old path failed in 10,498 of 20,000 cases, `mrev_edge` in none. Its worst arrival error, 1.5e-4 at `m = 30,000` and `t/t_min ~ 1e4`, is again at the one-ulp floor.
- **Everything else.**
  - Truth-set summaries are identical to the previous version.
  - Stress sets show 0 failures: 40,000 zero-revolution, 100,000 hyperbolic and 60,000 near-parabolic.
  - Multi-revolution 24,000 and 120,000: results unchanged, with status 2 only within 2.6e-9 of `t_min`.
  - Huge times: all 18,000 solved.
  - CTest passes 3/3 here and 9/9 in Lamberts-Comparison. On Linux x86-64 (Docker `gcc:14`, amd64), `test_excali` passes and a 10,000-case `abs(tau) > 0.9` stress run matches macOS.
  - No case with `abs(tau) <= 0.9` reaches `mrev_edge` in the 120,000-case set.
  - Rerun with both directions: the earlier 120,000-case multi-revolution and huge-time sets used direction +1 only, so the table path's long way had not been stress-tested. Results: 0 failures, status 2 only within 2.5e-9 of `t_min`, worst arrival error 1.9e-9 (2.8e-9 for the reference version); all 18,000 huge-time cases solved.
- **Speed.** `abs(tau) > 0.9` single solves take 2.22 µs, against 2.23 µs for the old path, which fails on 7% of them. A 32-node first ratio cost 2.93 µs: `libm` calls dominate each node, about 55 ns. Batch timing on the truth sets is unchanged: multirev 216.5 against 217.7 ns, zero_rev 182 to 183 ns.
- **Size.** `src/excali_solver.f90` is 1,675 lines (1,658 before).

## Single-precision tables; no iterated step in any starter (2026-11-06)

- **Tables in single precision.** Starters need about 0.1%; float32 rounding is 6e-8. Each value is written as the shortest decimal that rounds back to the same float32. The table source went from 606 KB to 231 KB, and the compiled tables from 163 KB to 81 KB. `mr_c` keeps its 1e-6 feasibility margin (interpolation 1.0e-7 plus rounding 6e-8). Truth and stress sets are unchanged. The only differences are status-2 flips at `t/t_min - 1 = 1e-9` and `1e-10`, the flag's threshold: 3 more of 780 adversarial near-minimum cases, flipping in both directions.
- **Body limits by table.** With `mu = 1/M` (short way) or `-1/M` (long way), the body limit is `Q(x) = (1 + mu (x + sin x))/cos^3(x/2) = 4 theta/M`. One family covers both ways.
  - Long-way minimum: `x_b = (8a/3)(1 + a^2 poly(a^2))`, `a = 1/M`, degree 4, within 2e-8. `Q(x_b) - 1` is evaluated directly, without cancellation; at a stationary point the error in `x_b` enters only to second order. Deriving it from `Tb_min = 1/(3 sin(x_b/2))` instead is first-order in that error, and at `m = 20000` it swamped the plateau depth (`mu^2` = 2.5e-10).
  - Inverse tables (8 `mu` points by 16), normalized so each value tends to a finite limit at both ends:
    - `x/xq` (short way) and `(x - x_b)/(xq - x_b)` (long way) up to `pi/2`, with `xq` the quadratic root about the origin or the minimum;
    - `cos(x/2) Q^(1/3)` beyond `pi/2`, shared by both ways in signed `mu`.

    Within 1.3e-4 (short way) and 2.7e-4 (long way) of `min(x - x_b, pi - x)`; 384 values.
  - Anchored body (long way, crowded side): the quadratic through `(0, A)` and the minimum, `x = x_b w/(1 + sqrt(1 - w))`. Within 1.3% of `x` (worst at `m = 1`). Its slope to the plateau is then explicit, so no inverse is needed.
  - Kappa-body (long way, far side): the layer term `cc/d^2` is bracketed by two independent inverses, without the layer and with it at `d*`. One linear interpolation of the residual gives the root within 0.037 of the bound. The previous version took three sequential inverses (0.007); two sequential inverses give 0.051, one gives 0.42.
- **Fold model by Ferrari.** `g(x) = (x-1)^2 (a2 + 2c/x + c/x^2)` exactly, so in `w = 1/x` the quartic `c w^4 + (a2 - 3c - dh) w^2 + 2(c - a2) w + a2 = 0` is already depressed. Ferrari's resolvent cubic (`depressed_max`) puts both roots on `w^2 - s w + k = 0`: the left root is `1/w+`, the right root `w+/k`. This matched bisection to 1e-10 (left) and 6e-6 (right) in 20,000 cases.
- **Starter error per kind** (Python mirror, 3,000 cases, error relative to the bound): unchanged except the kappa-body (0.007 to 0.037). Body 0.049, g-model 0.002, short-way fold cubic 0.124, all as before.
- **Results.**
  - Truth sets identical to the single-precision-table version.
  - zs2, hs, ps, ms (24k) and ms4d (120k): identical.
  - `abs(tau) > 0.9` stress, 40,000 cases, and the extreme variant: 0 failures. 37 and 6 cases are 10 times worse than the old path, against 38 and 7 before.
  - `mrev_edge` timing 2.35 to 2.16 us per solve, because the lookups are cheaper than the Newton steps.
  - Solver 1,669 lines (1,675 before).
- **Smaller multi-revolution root table.** `mr_q` moved from a uniform 32 `tau` × 16 `abs(sigma)` grid to 16 × 12: `tau = 0.9 sin(1.45 x)/sin(1.45)` with `x` uniform (one `asin` per lookup), and points uniform in `sqrt(abs(sigma))`. That is 4,224 values instead of 11,264, and more accurate. In 4,000 random cases (`m` from 1 to 1,000, `abs(sigma)` from `1e-2` to 1) the worst starter error relative to the bound is 1.4%, against 2.8% for the old grid; on another seed 0.97% against 1.1%. Both worst cases sit near `tau = -0.89`, `m = 2`. Uniform 16 × 16 reaches 5.4%, and uniform 16 × 8 7.4%.
  - Results: truth sets unchanged (adversarial category 17 maximum 2.3e-14 to 5.4e-14). 7,887 of the 10,000 multi-revolution truth cases change in the last bits, at most 2.4e-11 relative, with every summary the same. zs2, hs, ps, ms, ms4d and both `abs(tau) > 0.9` stress sets are identical. Timing unchanged.
- **Table totals.** 20,352 values (606 KB source, 163 KB compiled) became 13,696 values (155 KB source, 55 KB compiled).
- **Linux x86-64** (Docker gcc:14, amd64): CTest 3/3; the comparison repository's EXCALI test agrees with the contour solver in all five classes (0 bad); the `abs(tau) > 0.9` stress set (10,000 cases) and the extreme variant have no `-4`.

## One elliptic table (2026-11-06)

The zero-revolution elliptic starters `tab` (`abs(tau) <= 0.9`, split at `s = 1/2`) and `tabp` (`tau >= 0.75`, split at `v = pi/2`) became one table, `tabe`. It uses `tabp`'s layout for every `tau` in `[-0.9, 1-1e-6]`, on a grid uniform in `sign(tau)(1 - sqrt(1 - abs(tau)))`, so it is fine toward both ends. The pole side stores `(pi-v)/(s1 (K/T_half)^(1/3))`, `s1 = (T_half/t)^(1/3)`; the parabola side stores `v^2/(s2 (s2 + 2 e0))`, `s2 = (t-t_p)/(T_half-t_p)`. With the old one-sided grid (`1 - sqrt(1 - tau)`) the same layout reached 27% at `tau = -0.7`; the symmetric grid fixes that.

- **Size.** 40 `tau` points by 32: 2,560 values instead of 4,096. Ellipse starter formulas 8 to 6 (series, the table's two sides, `long_d`'s three pieces). Table totals 13,696 to 12,160 values.
- **Starter error** (Python mirror, 40,000 cases, relative to `min(v, pi-v)`):

  | `tau` | max error | rectangle |
  |---|---|---|
  | `[-0.9, -0.8)` | 5.9% | 8 nodes, bound/4 |
  | `[-0.8, -0.6)` | 1.5% | 8 nodes, bound/8 |
  | `[-0.6, 0.4)` | 0.14% | 3 nodes, bound/256 |
  | `[0.4, 1)` | 3.2e-5 | 3 nodes, bound/256 |

  The old `tab` was within 0.11% for `abs(tau) < 0.75` and 3.9% to 0.9, with 8 nodes at bound/8 from `-0.75`. `-0.75 <= tau < -0.6` now takes 8 nodes, and `-0.9 <= tau < -0.8` bound/4.
- **Results.** Truth sets: zero_rev, harness_mixed, multirev and infeasible identical; adversarial category 16 maximum 6.8e-16 to 9.3e-16, dev category 1 7.22e-15 to 7.24e-15. zs2 (40k), hs (100k), ps (60k): 0 failures, no velocity differs from the previous build by more than 1e-10. ms, ms4d, `abs(tau) > 0.9` stress and the extreme variant: identical.
- **Timing** (20,000 random zero-revolution ellipses, interleaved, machine under load): single 178-184 ns before, 177-179 after; batch 110-115 before, 111-117 after. Unchanged within noise.
- **Linux x86-64** (Docker gcc:14, amd64): the comparison repository's EXCALI test, 0 bad in all five classes; the `abs(tau) > 0.9` stress set and the extreme variant have no `-4`.
- **Alternatives studied** (Python mirror of the rectangle, same error measure):
  - Two rectangles with a 12×12×2 table (288 values) reach rounding level, but a second rectangle costs 41% single and 55% batch.
  - A fixed-node rectangle over all of `[0, pi]` in `X = v^2`, then a local rectangle, reaches rounding level for `abs(tau) <= 0.9` with no table, but fails toward `abs(tau) = 1`, where a singularity approaches the interval's end.
  - Izzo's starter and the bare pole and parabola asymptotes miss by about 100% of the bound near the ends, even with three rectangles.
  - No table reaches `tau < -0.9`: the long-way layer of width `eta^(2/3)` needs `long_d`.
