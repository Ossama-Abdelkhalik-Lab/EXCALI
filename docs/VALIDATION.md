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

The finite-rectangle extension uses a new implementation with adaptive coordinates and safeguarded corrections; it does not yet retain the original's lookup-table/vectorization optimizations. It is slower than ivLam on the local mixed benchmark. Original and extended timings must be reported separately.

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
