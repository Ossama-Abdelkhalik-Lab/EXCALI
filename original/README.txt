EXACT v2.2 vs EXACT v1 vs ivLam  (elliptic single-revolution Lambert problems)
==========================================================================

FILES
  exact_ser.F90, exact_vect.F90   Exact v2.2, serial / vectorized (= contour_study/gauss_legendre_*.F90)
  exact_v1.F90                    Exact v1 (code as received 2026-10-19, header comment added only)
  ivLam/                          ivLam runtime (ivLamRuntimeV2p50_739981p74401.f90, unchanged) + coefficient file
  benchmark_lambert.f90           speed + error histograms (both solvers called on the same data)
  tools/accuracy_compare.f90      |dV1| percentiles (50/90/99/99.9/99.99%)
  tools/vector_size_sweep.f90     time per solution vs problems per call (Q = 1 ... 5M)
  tools/v1_wrapper.f90            exposes Exact v1 to the drivers under the name rectsolution
  tools/plot_results.py           results/error_histogram.png, results/vector_size.png
  tools/generate_lambert.f90      test-data generator
  build.sh, run.sh                build and run everything (bash build.sh; bash run.sh; python tools/plot_results.py)
  results/data/                   raw outputs, prefix v22_ or v1_
Test data: ../fortranInterface/datalambert_fortran.bin (not copied, 1.1 GB). Single thread (OMP_NUM_THREADS=1).
Builds: serial  -O3 -fno-tree-vectorize -fno-tree-slp-vectorize
        vector  -O3 -march=native -fno-trapping-math

METHODS
  Exact v2.2   1. guess z from a 64 x 32 lookup table (Taylor step near z = 0)
               2. one precomputed thin rectangle (height 1e-7) around the guess, upper half only;
                  8 Gauss-Legendre points on the top edge + 1 per side; root = ratio of two weighted sums
               3. one Newton correction; retry on the correct rectangle if the guess picked the wrong one
               4. refuses the same 1,606 problems as ivLam (dt below the parabolic time); 180 deg fix for A
  Exact v1     1. coarse half circle, center 2 pi^2, radius 19.5 (covers z = 0.24 ... 39.24), 17 points
               2. fine half circle around the coarse root, radius <= 2, 9 points, complex sin/sinh/cos/cosh
               3. one complex-step Newton correction
               no lookup table, no retry, no check that the root lies inside the contour,
               OpenMP parallel loop (run here with 1 thread)
  ivLam        Russell's interpolation tree + corrections (reference)

ACCURACY  (|dV1|, 4,988,570 problems that ivLam solves, info = 0)
                         99.9%       (median)
  serial build
    Exact v2.2         3.75e-13    (2.99e-16)
    Exact v1           1.69e-01    (5.58e-16)
    ivLam              6.49e-13    (3.31e-16)
  vectorized build
    Exact v2.2         3.72e-13    (2.83e-16)
    Exact v1           1.69e-01    (5.14e-16)
    ivLam              6.31e-13    (3.16e-16)
  Exact v1: 1% of problems have error above 1.4e-5, 0.1% above 0.17.
  Refused problems: Exact v2.2 1,606 (same set as ivLam); Exact v1 0 (returns values for the 1,606
  problems ivLam refuses, errors 0.15 to Infinity).

SPEED  (5M problems in one call, ns per problem, median of 3 runs; ivLam timed in the same runs)
                         serial     vectorized
  Exact v2.2              98.8        42.6
  Exact v1               918.3       871.4
  ivLam                  ~140        ~145
  Exact v2.2 vs ivLam:   1.5x faster serial, 3.5x faster vectorized
  Exact v1 vs ivLam:     6.5x slower serial, 6.0x slower vectorized
  Exact v2.2 vs Exact v1: 9.3x faster serial, 20x faster vectorized

VECTOR SIZE  (ns per problem, problems split into calls of Q)
     Q      v2.2 ser  v2.2 vect   v1 ser   v1 vect   ivLam ser  ivLam vect
     1       130       139       5777     5660       139        139
     10       91        52       1375     1341       134        137
     100      95        38        952      964       144        137
     1000     89        35        897      870       138        134
     1M       86        38        892      890       136        138
  Exact v1 at small Q: about 5 us fixed cost per call, likely the OpenMP parallel region entered on every call.
