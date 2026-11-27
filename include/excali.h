/* EXCALI C interface (src/excali_c.f90). Started from gfortran -fc-prototypes,
 * then edited by hand; keep the two in step.
 *
 * Arrays of 3-vectors are double[n][3] (row i = case i), the same bytes as the
 * Fortran (3,n) arrays and as a C-contiguous NumPy array of shape (n,3): no
 * copy is made. Outputs are always written; check status before using them:
 *   0 solved; 1 coalesced multi-revolution roots at the minimum time;
 *   2 solved, low confidence; 5 time below the multi-revolution minimum;
 *  -1 invalid input; -2 degenerate geometry or missing antipodal plane;
 *  -3 nonunique collinear geometry; -4 numerical failure.
 * Failed cases return NaN velocities. direction = +1 follows r1 x r2, -1 the
 * opposite; nrev >= 0; branch = -1 / +1 the lower / higher multi-revolution
 * root (either value for nrev = 0). */
#ifndef EXCALI_H
#define EXCALI_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* One case in any consistent units: mu in distance^3/time^2. */
void excali_solve(const double r1[3], const double r2[3], double tof, double mu,
                  int direction, int nrev, int branch,
                  double v1[3], double v2[3], int *status);

/* As excali_solve, with the transfer-plane normal for exactly antipodal
 * endpoints (perpendicular to both position vectors). */
void excali_solve_plane(const double r1[3], const double r2[3], double tof, double mu,
                        int direction, int nrev, int branch, const double plane_normal[3],
                        double v1[3], double v2[3], int *status);

/* n cases with mu = 1 (the Fortran excali_batch). half[i] is the geometry
 * flag: -1 degenerate, 0 ordinary, 1 near-collinear. */
void excali_batch(int n, const double (*r1)[3], const double (*r2)[3], const double *tof,
                  const int *direction, const int *nrev, const int *branch,
                  double (*v1)[3], double (*v2)[3], int *status, int *half);

/* n cases sharing one mu. Time is rescaled by sqrt(mu) in a temporary copy of
 * tof; results agree with per-case excali_solve to rounding (statuses
 * identical). Extreme units that overflow tof*sqrt(mu) return status -1:
 * use excali_solve there. */
void excali_batch_mu(int n, const double (*r1)[3], const double (*r2)[3], const double *tof,
                     double mu, const int *direction, const int *nrev, const int *branch,
                     double (*v1)[3], double (*v2)[3], int *status, int *half);

/* Diagnostic: the address Fortran receives for x (tests use it to check
 * that arrays are passed without a copy). */
intptr_t excali_c_data_address(const double *x);

#ifdef __cplusplus
}
#endif

#endif
