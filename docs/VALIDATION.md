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
