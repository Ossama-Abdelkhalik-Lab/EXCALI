/* The C interface (include/excali.h) against direct Fortran calls of the
 * same routines (tests/excali_c_ref.f90): identical statuses and bit-identical
 * velocities for single cases and batches, on random zero- and
 * multi-revolution transfers, extreme units, antipodal endpoints with a plane,
 * and invalid inputs. Also: excali_batch_mu against per-case excali_solve
 * (statuses identical, velocities to rounding) and the no-copy address check. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "excali.h"

void ref_solve(const double *r1, const double *r2, double tof, double mu, int direction,
               int nrev, int branch, double *v1, double *v2, int *status);
void ref_solve_plane(const double *r1, const double *r2, double tof, double mu, int direction,
                     int nrev, int branch, const double *pn, double *v1, double *v2, int *status);
void ref_batch(int n, const double *r1, const double *r2, const double *tof, const int *direction,
               const int *nrev, const int *branch, double *v1, double *v2, int *status, int *half);

static uint64_t rng = 0x9E3779B97F4A7C15ull;
static double uni(void) { /* xorshift64*, uniform in [0,1) */
  rng ^= rng >> 12; rng ^= rng << 25; rng ^= rng >> 27;
  return (double)((rng * 2685821657736338717ull) >> 11) * 0x1.0p-53;
}

/* Bitwise equality, with any NaN equal to any NaN. */
static int same(double a, double b) {
  if (isnan(a) && isnan(b)) return 1;
  return memcmp(&a, &b, sizeof a) == 0;
}

#define NC 600
static double R1[NC][3], R2[NC][3], TOF[NC], MU[NC];
static int DIR[NC], NREV[NC], BR[NC];

static void make_cases(void) {
  const double pi = acos(-1.0);
  for (int i = 0; i < NC; i++) {
    double r1 = 0.5 + 2.5 * uni(), r2 = 0.5 + 2.5 * uni();
    double th = 0.05 + (2 * pi - 0.1) * uni(), ph = 2 * pi * uni();
    double e1[3] = {cos(ph), sin(ph), 0}, e2[3] = {-sin(ph) * cos(1.0), cos(ph) * cos(1.0), sin(1.0)};
    for (int k = 0; k < 3; k++) {
      R1[i][k] = r1 * e1[k];
      R2[i][k] = r2 * (cos(th) * e1[k] + sin(th) * e2[k]);
    }
    DIR[i] = uni() < 0.5 ? 1 : -1;
    MU[i] = 0.25 + 4 * uni();
    double s = 0.5 * (r1 + r2 + sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * cos(th)));
    double unit = sqrt(s * s * s / MU[i]); /* time scale */
    if (i < 300) {
      NREV[i] = 0; BR[i] = -1;
      TOF[i] = unit * pow(10, -1.5 + 2.5 * uni()); /* hyperbolic to long ellipses */
    } else if (i < 520) {
      NREV[i] = 1 + (int)(5 * uni()); BR[i] = uni() < 0.5 ? 1 : -1;
      TOF[i] = unit * NREV[i] * pi / sqrt(2.0) * (0.9 + 3 * uni()); /* some infeasible */
    } else {
      NREV[i] = 0; BR[i] = 1; TOF[i] = unit; /* overwritten below: special cases */
    }
  }
  /* Invalid and special inputs (indices 520..599). */
  int j = 520;
  TOF[j++] = -1;                       /* tof <= 0 */
  TOF[j++] = 0;
  TOF[j] = NAN; j++;
  R1[j][0] = NAN; j++;                 /* non-finite position */
  MU[j++] = 0;                         /* mu <= 0 */
  MU[j++] = -2;
  MU[j] = INFINITY; j++;
  DIR[j++] = 0;                        /* bad direction */
  NREV[j++] = -1;                      /* bad revolutions */
  BR[j] = 0; NREV[j] = 2; j++;         /* bad branch */
  for (int k = 0; k < 3; k++) R1[j][k] = 0; j++;       /* zero radius */
  for (int k = 0; k < 3; k++) R2[j][k] = 2 * R1[j][k]; j++;   /* same ray */
  for (int k = 0; k < 3; k++) R2[j][k] = -R1[j][k]; j++;      /* antipodal, no plane */
  for (; j < 560; j++) {               /* extreme units, mu*tof^2 spanning 1e+-200 */
    double sc = pow(10, -150 + 300 * uni());
    for (int k = 0; k < 3; k++) { R1[j][k] *= sc; R2[j][k] *= sc; }
    MU[j] *= pow(10, -100 + 200 * uni());
    double s0 = sqrt(R1[j][0] * R1[j][0] + R1[j][1] * R1[j][1] + R1[j][2] * R1[j][2]);
    TOF[j] = (0.5 + 3 * uni()) * sqrt(s0 * s0 * s0 / MU[j]);
  }
  for (; j < NC; j++) {                /* near-parabolic and near-minimum-time */
    if (j % 2) { NREV[j] = 1 + j % 4; BR[j] = j % 3 ? 1 : -1; }
  }
}

int main(void) {
  int bad = 0;
  make_cases();

  /* Single cases: C entry vs direct Fortran call. */
  int nstat[8] = {0};
  for (int i = 0; i < NC; i++) {
    double a1[3], a2[3], b1[3], b2[3];
    int sa, sb;
    excali_solve(R1[i], R2[i], TOF[i], MU[i], DIR[i], NREV[i], BR[i], a1, a2, &sa);
    ref_solve(R1[i], R2[i], TOF[i], MU[i], DIR[i], NREV[i], BR[i], b1, b2, &sb);
    int ok = sa == sb;
    for (int k = 0; k < 3; k++) ok = ok && same(a1[k], b1[k]) && same(a2[k], b2[k]);
    if (!ok) { bad++; printf("solve mismatch case %d: status %d vs %d\n", i, sa, sb); }
    if (sa >= -4 && sa <= 2) nstat[sa + 4]++; else if (sa == 5) nstat[7]++;
  }
  printf("single cases %d, statuses -4..2,5: %d %d %d %d %d %d %d %d\n", NC,
         nstat[0], nstat[1], nstat[2], nstat[3], nstat[4], nstat[5], nstat[6], nstat[7]);
  if (nstat[4] < 400) { bad++; printf("too few solved cases (%d)\n", nstat[4]); }

  /* Antipodal endpoints with a transfer plane. */
  {
    double r1[3] = {1, 0, 0}, r2[3] = {-1.5, 0, 0}, pn[3] = {0, 0, 1}, a1[3], a2[3], b1[3], b2[3];
    int sa, sb, sn;
    excali_solve_plane(r1, r2, 4.0, 1.0, 1, 0, -1, pn, a1, a2, &sa);
    ref_solve_plane(r1, r2, 4.0, 1.0, 1, 0, -1, pn, b1, b2, &sb);
    excali_solve(r1, r2, 4.0, 1.0, 1, 0, -1, b1, b2, &sn);
    ref_solve_plane(r1, r2, 4.0, 1.0, 1, 0, -1, pn, b1, b2, &sb);
    int ok = sa == 0 && sb == 0 && sn == -2;
    for (int k = 0; k < 3; k++) ok = ok && same(a1[k], b1[k]) && same(a2[k], b2[k]);
    if (!ok) { bad++; printf("antipodal plane mismatch: %d %d %d\n", sa, sb, sn); }
  }

  /* Batches (mu = 1): C entry vs direct Fortran call, several sizes. */
  static double T1[NC], A1[NC][3], A2[NC][3], B1[NC][3], B2[NC][3];
  static int SA[NC], SB[NC], HA[NC], HB[NC];
  for (int i = 0; i < NC; i++) T1[i] = TOF[i] * sqrt(MU[i]);
  int sizes[4] = {0, 1, 32, NC};
  for (int q = 0; q < 4; q++) {
    int n = sizes[q];
    for (int i = 0; i < NC; i++) { SA[i] = SB[i] = 99; HA[i] = HB[i] = 99; }
    excali_batch(n, R1, R2, T1, DIR, NREV, BR, A1, A2, SA, HA);
    if (n > 0) ref_batch(n, &R1[0][0], &R2[0][0], T1, DIR, NREV, BR, &B1[0][0], &B2[0][0], SB, HB);
    for (int i = 0; i < NC; i++) {
      int ok = i < n ? (SA[i] == SB[i] && HA[i] == HB[i]) : (SA[i] == 99 && HA[i] == 99);
      for (int k = 0; i < n && k < 3; k++) ok = ok && same(A1[i][k], B1[i][k]) && same(A2[i][k], B2[i][k]);
      if (!ok) { bad++; printf("batch n=%d mismatch case %d: %d vs %d\n", n, i, SA[i], SB[i]); break; }
    }
  }

  /* excali_batch_mu against the mu = 1 batch scaled here (bit-identical:
   * the same two multiplications) and against per-case excali_solve. */
  {
    double mu = 398600.4418, sm = sqrt(mu), tm[NC];
    for (int i = 0; i < NC; i++) tm[i] = T1[i] / sm;
    excali_batch_mu(NC, R1, R2, tm, mu, DIR, NREV, BR, A1, A2, SA, HA);
    for (int i = 0; i < NC; i++) T1[i] = tm[i] * sm;
    ref_batch(NC, &R1[0][0], &R2[0][0], T1, DIR, NREV, BR, &B1[0][0], &B2[0][0], SB, HB);
    double worst = 0;
    for (int i = 0; i < NC; i++) {
      int ok = SA[i] == SB[i];
      for (int k = 0; k < 3; k++) {
        double b1 = SB[i] >= 0 && SB[i] <= 2 ? B1[i][k] * sm : B1[i][k];
        double b2 = SB[i] >= 0 && SB[i] <= 2 ? B2[i][k] * sm : B2[i][k];
        ok = ok && same(A1[i][k], b1) && same(A2[i][k], b2);
      }
      if (!ok) { bad++; printf("batch_mu scaling mismatch case %d\n", i); break; }
      double c1[3], c2[3];
      int sc;
      excali_solve(R1[i], R2[i], tm[i], mu, DIR[i], NREV[i], BR[i], c1, c2, &sc);
      if (sc != SA[i]) { bad++; printf("batch_mu vs solve status case %d: %d vs %d\n", i, SA[i], sc); continue; }
      if (sc >= 0 && sc <= 2) {
        double d = 0, nv = 0;
        for (int k = 0; k < 3; k++) {
          d = fmax(d, fmax(fabs(A1[i][k] - c1[k]), fabs(A2[i][k] - c2[k])));
          nv = fmax(nv, fabs(c1[k]));
        }
        if (d / nv > worst) worst = d / nv;
      }
    }
    printf("batch_mu vs per-case solve: worst relative velocity difference %.2e\n", worst);
    if (worst > 1e-9) { bad++; printf("batch_mu disagrees with excali_solve beyond rounding\n"); }
    /* Invalid mu: every case -1 with NaN outputs. */
    excali_batch_mu(3, R1, R2, tm, -1.0, DIR, NREV, BR, A1, A2, SA, HA);
    for (int i = 0; i < 3; i++)
      if (SA[i] != -1 || !isnan(A1[i][0])) { bad++; printf("invalid mu not rejected\n"); break; }
  }

  /* No copy: Fortran sees the caller's addresses. */
  if (excali_c_data_address(&R1[5][0]) != (intptr_t)&R1[5][0]) { bad++; printf("address differs\n"); }

  printf("%s (%d failures)\n", bad ? "FAIL" : "PASS", bad);
  return bad ? 1 : 0;
}
