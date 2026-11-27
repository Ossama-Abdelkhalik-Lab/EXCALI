// EXCALI zero-revolution path as __host__ __device__-style C++ (prototype).
//
// A line-by-line port of excali_one for nrev = 0 from src/excali_solver.f90 and
// src/excali_stages.inc (main f5ca647): geometry, zero_rev_start (parabola,
// near-parabolic series, tabe / long_d ellipses, tabh and large-w hyperbolas
// including the xi coordinate), fast_coordinate, rect_ratio (kinds 0, 1, 2, -1,
// -2, with the exact node-hit guard), fast_accept and velocities. Multiple
// revolutions are not ported (status -6).
//
// The same source compiles as plain C++ (host test against Fortran) and as
// CUDA device code (EXC_FN = __device__ under nvcc). Tables come from
// excali_tables_cuda.h, generated bit-exact from the Fortran constants; the
// trig tables are the symmetric positive half.
#pragma once
#include <cmath>
#include <cstdint>
#include <cstring>
#include <climits>

#ifdef __CUDACC__
#define EXC_FN __device__ __forceinline__
#ifndef EXCALI_TABLE
#ifdef EXCALI_TABLES_GLOBAL
#define EXCALI_TABLE __device__ const        // global memory (read-only data cache)
#else
#define EXCALI_TABLE __constant__ const      // constant memory
#endif
#endif
#else
#define EXC_FN static inline
#endif
#include "excali_tables_cuda.h"

namespace exc {

constexpr double pi = EXC_PI, sqrt2 = EXC_SQRT2, rsqrt2 = EXC_RSQRT2;
constexpr double eps = 2.220446049250313e-16;          // epsilon(1._dp)
constexpr double huge_ = 1.7976931348623157e308;       // huge(1._dp)
constexpr double tiny_ = 2.2250738585072014e-308;      // tiny(1._dp)
constexpr int kmax = 62, np_ser = 18, tabe_nt = 40, tabe_ns = 32, tab_n = 32, bb_n = 6;

EXC_FN int64_t dbits(double x) {
#ifdef __CUDA_ARCH__
    return __double_as_longlong(x);
#else
    int64_t i; std::memcpy(&i, &x, 8); return i;
#endif
}
EXC_FN double bitsd(int64_t i) {
#ifdef __CUDA_ARCH__
    return __longlong_as_double(i);
#else
    double x; std::memcpy(&x, &i, 8); return x;
#endif
}
EXC_FN double sgn(double a, double b) { return std::signbit(b) ? -std::fabs(a) : std::fabs(a); }
EXC_FN bool finite3(const double v[3]) { return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]); }
EXC_FN double dot3(const double a[3], const double b[3]) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
EXC_FN double norm3(const double a[3]) { return std::sqrt(dot3(a, a)); }
EXC_FN void cross3(const double u[3], const double v[3], double w[3]) {
    w[0] = u[1] * v[2] - u[2] * v[1]; w[1] = u[2] * v[0] - u[0] * v[2]; w[2] = u[0] * v[1] - u[1] * v[0];
}

// Half-table node layout: node j (0-based, 0..9) -> positive index and sign.
EXC_FN int hidx(int j) { return j < 4 ? j : (j < 8 ? j - 4 : 4); }
EXC_FN double hsgn(int j) { return (j < 4 || j == 8) ? 1.0 : -1.0; }

struct Frame {
    double scale, rn1, rn2, e1[3], e2[3], t1v[3], t2v[3], cs, ss, s1, s2, p, q, a, b, y0, yend, ypi, t;
};

EXC_FN double pow2m(int k) { return bitsd(((int64_t)(1023 - k)) << 52); }
EXC_FN int kexp(double x) {
    if (x >= tiny_) return 1023 - (int)((uint64_t)dbits(x) >> 52);
    if (x > 0) return -std::ilogb(x);                   // 1 - exponent(x)
    return INT_MAX;
}
EXC_FN double rlog(double x) {
    const double ln2 = 0.6931471805599453;               // log(2._dp)
    if (!(x >= tiny_ && x <= huge_)) return std::log(x);
    int64_t i = dbits(x), e = (i >> 52) - 1023;
    double m = bitsd((i & 4503599627370495LL) | 4607182418800017408LL);
    if (m > sqrt2) { m = m / 2; e = e + 1; }
    double z = (m - 1) / (m + 1), q = z * z;
    return e * ln2 + 2 * z * (1 + q * (1.0 / 3 + q * (.2 + q * (1.0 / 7 + q * (1.0 / 9)))));
}
EXC_FN double rcbrt(double x) {                          // coarse (one Halley step)
    if (!(x >= tiny_ && x <= huge_)) return std::cbrt(x);
    double y = bitsd(dbits(x) / 3 + 3071306045661749710LL);
    double w = 1 / y, r = (x * w) * (w * w);
    return y * (1 + 2 * r) / (2 + r);
}
EXC_FN void lagrange4(double p, double w[4]) {
    double p0 = p, p1 = p - 1, p2 = p - 2, p3 = p - 3;
    w[0] = -(1.0 / 6) * (p1 * p2) * p3; w[1] = .5 * (p0 * p2) * p3;
    w[2] = -.5 * (p0 * p1) * p3;       w[3] = (1.0 / 6) * (p0 * p1) * p2;
}
EXC_FN void shch(double x, double &s, double &c) {
    if (std::fabs(x) < .5) {
        double q = x * x; const double *rf = exc_rf;
        s = x * (1 + q * (rf[2] + q * (rf[4] + q * (rf[6] + q * (rf[8] + q * (rf[10] + q * rf[12]))))));
        c = 1 + q * (rf[1] + q * (rf[3] + q * (rf[5] + q * (rf[7] + q * (rf[9] + q * (rf[11] + q * rf[13]))))));
    } else if (std::fabs(x) < 700) {
        double e = std::exp(x), q = 1 / e; s = (e - q) / 2; c = (e + q) / 2;
    } else { s = std::sinh(x); c = std::cosh(x); }
}
EXC_FN void sincos_half(double c, double &s, double &co) {
    const double pio2_hi = 1.5707963267948966, pio2_lo = 6.123233995736766e-17;
    if (!(c >= 0 && c <= pi)) { s = std::sin(c / 2); co = std::cos(c / 2); return; }
    const double *rf = exc_rf;
    double y = c / 2;
    if (y > pio2_hi / 2) y = (pio2_hi - y) + pio2_lo;
    double z = y * y;
    double ps = y + y * z * (-rf[2] + z * (rf[4] + z * (-rf[6] + z * (rf[8] + z * (-rf[10] + z * (rf[12] - z * rf[14]))))));
    double pc = 1 + z * (-rf[1] + z * (rf[3] + z * (-rf[5] + z * (rf[7] + z * (-rf[9] + z * (rf[11] + z * (-rf[13] + z * rf[15])))))));
    if (c / 2 > pio2_hi / 2) { s = pc; co = ps; } else { s = ps; co = pc; }
}
EXC_FN void rotate(int kind1, double x, double h, double &hs, double &hc) {
    if (kind1 == 0) return;
    double s, c, q, r; const double *rf = exc_rf;
    if (std::fabs(h) < .25) {
        q = h * h; if (kind1 > 0) q = -q;
        s = h * (1 + q * (rf[2] + q * (rf[4] + q * (rf[6] + q * (rf[8] + q * (rf[10] + q * rf[12]))))));
        c = 1 + q * (rf[1] + q * (rf[3] + q * (rf[5] + q * (rf[7] + q * (rf[9] + q * (rf[11] + q * rf[13]))))));
    } else if (kind1 > 0) { s = std::sin(h); c = std::cos(h); }
    else shch(h, s, c);
    q = hs; r = hc;
    hs = q * c + hc * s;
    if (kind1 > 0) hc = hc * c - q * s; else hc = hc * c + q * s;
    if (std::fabs(hs) < std::fabs(q) / 4 || std::fabs(hc) < std::fabs(r) / 4) {
        if (kind1 > 0) { hs = std::sin(x / 2); hc = std::cos(x / 2); } else shch(x / 2, hs, hc);
    }
}
EXC_FN double depressed_max(double p, double q) {
    double m = std::sqrt(-p / 3), c = -q / (2 * (m * m * m));
    if (std::fabs(c) <= 1) return 2 * m * std::cos(std::acos(c) / 3);
    return sgn(2 * m * std::cosh(std::log(std::fabs(c) + std::sqrt(c * c - 1)) / 3), c);
}
EXC_FN double cubic_pos(double e, double f) {
    double s = std::cbrt(f), chi = e / s;
    if (chi > 0) return s / depressed_max(-chi, -1.0);
    if (chi < 0) return s * (depressed_max(-chi * chi / 3, 2 * (chi * chi * chi) / 27 - 1) - chi / 3);
    return s;
}
EXC_FN void fhyp_d(double w, double sh, double ch, double &f, double &f1) {
    if (w < EXC_XSER) {
        double xx = -w * w; f = exc_pser2[0]; f1 = 0;
        for (int j = 1; j < np_ser; ++j) { f1 = f1 * xx + f; f = f * xx + exc_pser2[j]; }
        f1 = -2 * w * f1;
    } else {
        f = (sh * ch - w) * rsqrt2 / (sh * sh * sh); f1 = (sqrt2 - 3 * ch * f) / sh;
    }
}
EXC_FN void near_parabola_x(double y0, double a, double b, double t, double tp, double &xx, double &rc) {
    double g0 = exc_pser2[np_ser - 1], g1 = exc_pser2[np_ser - 2], g2 = exc_pser2[np_ser - 3];
    double q0 = std::sqrt(y0), q1 = b / (4 * q0), q2 = (-b / 24 - q1 * q1) / (2 * q0);
    double p0 = y0 * g0 + a, p1 = y0 * g1 + b / 2 * g0, p2 = y0 * g2 + b / 2 * g1 - b / 24 * g0;
    double t1 = q0 * p1 + q1 * p0, t2 = q0 * p2 + q1 * p1 + q2 * p0;
    double d = t - tp;
    rc = t2 * d / (t1 * t1);
    xx = d / t1 * (1 - rc);
}
EXC_FN float tabe_at(int i, int j, int side) { return exc_tabe[(i - 1) + tabe_nt * (j - 1) + tabe_nt * tabe_ns * (side - 1)]; }
EXC_FN float tabh_at(int i, int j, int h) { return exc_tabh[(i - 1) + tab_n * (j - 1) + tab_n * tab_n * (h - 1)]; }
EXC_FN double table_ve(double tau, double s, double sc, int side) {
    double u = sgn(1 - std::sqrt(1 - std::fabs(tau)), tau);
    if (u > EXC_TABE_Z1) u = EXC_TABE_Z1;
    if (!(u >= EXC_TABE_Z0 && s >= 0 && s <= 1)) return 0;
    double wt[4], ws[4], r[4];
    u = (u - EXC_TABE_Z0) * ((tabe_nt - 1) / (EXC_TABE_Z1 - EXC_TABE_Z0));
    int i = (int)u - 1; if (i < 0) i = 0; if (i > tabe_nt - 4) i = tabe_nt - 4;
    lagrange4(u - i, wt);
    u = s * (tabe_ns - 1);
    int j = (int)u - 1; if (j < 0) j = 0; if (j > tabe_ns - 4) j = tabe_ns - 4;
    lagrange4(u - j, ws);
    for (int c = 0; c < 4; ++c)
        r[c] = wt[0] * tabe_at(i + 1, j + 1 + c, side) + wt[1] * tabe_at(i + 2, j + 1 + c, side)
             + wt[2] * tabe_at(i + 3, j + 1 + c, side) + wt[3] * tabe_at(i + 4, j + 1 + c, side);
    u = (ws[0] * r[0] + ws[1] * r[1]) + (ws[2] * r[2] + ws[3] * r[3]);
    double v;
    if (side == 1) v = pi - s * sc * u;
    else { double z = u * sc; v = std::sqrt(z > 0 ? z : 0.0); }
    if (!(v > 0 && v < pi)) v = 0;
    return v;
}
EXC_FN double long_d(double th, double eta) {
    double dl = th - pi / 4, x, u, h, z, d;
    if (dl > 0) {
        x = 4 * dl / pi;
        if (x < 1e-3) x = x * (2.0 / 3 - x * (1.0 / 9 - x * (4.0 / 81)));
        else x = std::exp(std::log1p(x) * (2.0 / 3)) - 1;
        d = eta / std::sqrt(x);
        if (d < 2 * eta) return d;
    }
    d = cubic_pos(2 * dl, .75 * pi * eta * eta);
    if (d < .1) return d;
    u = (th - 1.0 / 3) / (pi / 4 - 1.0 / 3); if (u < 0) u = 0; if (u > 1) u = 1;
    const double *bb = exc_bb;
    if (u <= .5) {
        z = 4 * u - 1; h = bb[bb_n - 1];
        for (int j = bb_n - 2; j >= 0; --j) h = h * z + bb[j];
        double hu = h * u; d = pi - std::sqrt(hu > 0 ? hu : 0.0);
    } else {
        z = 4 * u - 3; h = bb[bb_n + bb_n - 1];
        for (int j = bb_n - 2; j >= 0; --j) h = h * z + bb[bb_n + j];
        d = h * (1 - u);
    }
    return d;
}
EXC_FN double table_w(double tau, double u, double wb) {
    if (!(u >= EXC_TABH_U0 && u <= .99)) return 0;
    int h = tau > 0 ? 1 : 2;
    double wt[4], ws[4], cl[4], s, c;
    double g = std::sqrt(std::fabs(tau)) * ((tab_n - 1) / EXC_TABH_S1);
    if (g > tab_n - 1.0) g = tab_n - 1.0;
    int i = (int)g - 1; if (i < 0) i = 0; if (i > tab_n - 4) i = tab_n - 4;
    lagrange4(g - i, wt);
    g = (u - EXC_TABH_U0) * ((tab_n - 1) / (1 - EXC_TABH_U0));
    int j = (int)g - 1; if (j < 0) j = 0; if (j > tab_n - 4) j = tab_n - 4;
    lagrange4(g - j, ws);
    for (int k = 0; k < 4; ++k)
        cl[k] = wt[0] * tabh_at(i + 1, j + 1 + k, h) + wt[1] * tabh_at(i + 2, j + 1 + k, h)
              + wt[2] * tabh_at(i + 3, j + 1 + k, h) + wt[3] * tabh_at(i + 4, j + 1 + k, h);
    double r = (ws[0] * cl[0] + ws[1] * cl[1]) + (ws[2] * cl[2] + ws[3] * cl[3]);
    double l2 = 1 - rlog(u); l2 = l2 * l2;
    if (wb > 0) {
        double z = wb * wb / l2;
        if (z < 1) { shch(z / 2, s, c); l2 = l2 * 2 * s * (c - s); }
        else l2 = l2 * (1 - std::exp(-z));
    }
    double rr = r * (1 - u) * l2;
    double w = std::sqrt(rr > 0 ? rr : 0.0);
    if (wb > 0 && !(w < wb)) w = 0;
    if (!(w > 0)) w = 0;
    return w;
}

// Transfer frame normalized by the larger radius; status -4 means usable.
EXC_FN void geometry(const double r1[3], const double r2[3], double tof, int direction, Frame &g,
                     int &status, int &half, const double *plane_normal) {
    double h[3], h2[3], hn, dr, d12[3], s12[3];
    half = 0; status = -1;
    g.rn1 = std::sqrt(dot3(r1, r1)); g.rn2 = std::sqrt(dot3(r2, r2));
    if (!(g.rn1 > 0 && g.rn2 > 0) || !std::isfinite(g.rn1 + g.rn2)) return;
    dr = g.rn1 - g.rn2;
    if (std::fabs(dr) < .5 * (g.rn1 > g.rn2 ? g.rn1 : g.rn2)) {
        for (int k = 0; k < 3; ++k) { d12[k] = r1[k] - r2[k]; s12[k] = r1[k] + r2[k]; }
        dr = dot3(d12, s12) / (g.rn1 + g.rn2);
    }
    for (int k = 0; k < 3; ++k) { g.e1[k] = r1[k] / g.rn1; g.e2[k] = r2[k] / g.rn2; }
    if (g.rn1 >= g.rn2) { g.scale = g.rn1; g.rn2 = g.rn2 / g.rn1; g.rn1 = 1; }
    else { g.scale = g.rn2; g.rn1 = g.rn1 / g.rn2; g.rn2 = 1; }
    cross3(g.e1, g.e2, h); hn = norm3(h);
    status = -2;
    if (hn < 32 * eps) {
        half = -1;
        if (!plane_normal) return;
        for (int k = 0; k < 3; ++k) h2[k] = g.e1[k] + g.e2[k];
        if (norm3(h2) > 32 * eps) { status = -3; return; }
        if (!finite3(plane_normal)) return;
        for (int k = 0; k < 3; ++k) h[k] = plane_normal[k];
        hn = norm3(h);
        if (hn <= 0 || !std::isfinite(hn)) return;
        double hh[3] = {h[0] / hn, h[1] / hn, h[2] / hn};
        double m1 = std::fabs(dot3(hh, g.e1)), m2 = std::fabs(dot3(hh, g.e2));
        if ((m1 > m2 ? m1 : m2) > 64 * eps) return;
        half = 1;
    }
    if (hn < 1e-4) half = 1;
    double f = direction / hn;
    for (int k = 0; k < 3; ++k) h[k] = h[k] * f;
    for (int k = 0; k < 3; ++k) h2[k] = g.e1[k] + g.e2[k];
    g.cs = std::sqrt(dot3(h2, h2)) / 2;
    for (int k = 0; k < 3; ++k) h2[k] = g.e2[k] - g.e1[k];
    g.ss = std::sqrt(dot3(h2, h2)) / 2;
    cross3(h, g.e1, g.t1v); cross3(h, g.e2, g.t2v);
    g.s1 = std::sqrt(g.rn1); g.s2 = std::sqrt(g.rn2);
    g.p = g.s1 * g.s2; g.a = direction * sqrt2 * g.p * g.cs; g.b = sqrt2 * g.a;
    g.q = (g.rn1 - g.rn2) / (g.s1 + g.s2);
    if ((g.rn1 < g.rn2 ? g.rn1 : g.rn2) > .5) g.q = (dr / g.scale) / (g.s1 + g.s2);
    g.yend = g.q * g.q + 2 * g.p * g.ss * g.ss / (1 + g.cs);
    g.y0 = g.yend; if (direction < 0) g.y0 = g.q * g.q + 2 * g.p * (1 + g.cs);
    g.ypi = g.yend; if (g.b >= 0) g.ypi = g.y0 + 2 * g.b;
    g.t = tof / (g.scale * std::sqrt(g.scale));
    status = -4;
}

EXC_FN void zero_rev_start(const Frame &g, int &k, double &x, double &wb, bool &fin, bool &nearp) {
    double tp = std::sqrt(g.y0) * (g.y0 * sqrt2 / 3 + g.a), u = 0, tau, dc, sl, sw, r0, r1, dt;
    k = 0; x = 0; wb = 0; fin = false; nearp = false;
    if (std::fabs(g.t - tp) <= 4 * eps * tp) return;
    nearp = std::fabs(g.t - tp) <= .031 * tp;
    k = 1;
    if (g.t < tp) {
        k = -1; u = g.t / tp;
        if (g.b > 0) { sw = std::sqrt(g.y0 / (2 * g.b)); wb = 2 * std::log1p(sw + sw * sw / (1 + std::sqrt(1 + sw * sw))); }
    }
    if (nearp) {
        near_parabola_x(g.y0, g.a, g.b, g.t, tp, u, dt);
        if (u * k > 0) x = std::sqrt(std::fabs(u));
        fin = std::fabs(dt) < 1e-8;
        return;
    }
    if (k < 0 && u < EXC_TABH_U0) {
        sw = g.y0 + g.b; double q = sw * sw * g.b;
        dc = std::sqrt(q * q + 8 * g.t * g.t * (sw * sw * sw));
        if (g.b >= 0) dc = 4 * (sw * sw * sw) / (sw * sw * g.b + dc);
        else dc = (dc - sw * sw * g.b) / (2 * g.t * g.t);
        x = std::log(dc);
        if (g.b > 0) {
            sl = g.t * dc / (sqrt2 * sw); sl = sl * sl;
            if (sl < g.y0 / 2) {
                r0 = std::sqrt(g.y0 * g.ypi); r1 = std::sqrt((g.y0 - sl) * (g.ypi - sl));
                k = -2; x = std::log1p((sl + sl * (2 * sw - sl) / (r0 + r1)) / (sw - sl + r1));
            } else { double c99 = .99 * wb; if (c99 < x) x = c99; }
        }
        return;
    }
    double s12 = g.s1 + g.s2;
    dc = std::sqrt(g.q * g.q * (s12 * s12) + 4 * g.p * g.p * g.ss * g.ss); sl = (g.rn1 + g.rn2 + dc) / 2;
    tau = g.a / (sqrt2 * sl) / (1 + std::sqrt(dc / sl));
    if (k < 0) { x = table_w(tau, u, wb); return; }
    if (tau < -.9) {
        x = long_d(g.t / ((-g.b) * std::sqrt(-g.b)), std::sqrt(2 * g.yend / (-g.b)));
        if (x > .25 * pi) x = pi - x; else k = 2;
        return;
    }
    sw = g.y0 + g.b;
    sw = std::sqrt(sw) * (sw * (pi / (2 * sqrt2)) + g.a);
    if (g.t >= sw) {
        u = rcbrt(sw / g.t); x = table_ve(tau, u, rcbrt(pi * g.ypi * std::sqrt(g.ypi) / (sqrt2 * sw)), 1);
    } else {
        u = (g.t - tp) / (sw - tp); x = table_ve(tau, u, u * (u + 2 * tp / (sw - tp)), 2);
    }
    if (x > .75 * pi) { x = pi - x; k = 2; }
}

EXC_FN void fast_coordinate(int k, double x, double y0, double b, double wb, int &kind1, double &c,
                            double &el, double &rb, double &hs, double &hc) {
    hs = 0; hc = 1; kind1 = k; c = x;
    double ac;
    if ((k == 1 || k == -1) && x < EXC_XSER) {
        kind1 = 0; c = k * x * x; ac = std::fabs(c);
        el = std::fmin(ac, EXC_XWIN - ac);
        if (k > 0) el = std::fmin(2 * ac, EXC_XWIN - ac);
        if (k > 0 && b > 0) el = std::fmin(el, c + 2 * y0 / (b + y0 / 2));
        if (k < 0 && b > 0) el = std::fmin(el, wb * wb - ac);
    } else if (k > 0) {
        el = std::fmin(c, pi - c); sincos_half(c, hs, hc);
    } else if (k == -2) {
        el = std::fmin(std::fmin(c, 2 * (wb - c)), pi); shch(c / 2, hs, hc);
    } else {
        el = std::fmin(c, pi);
        if (b > 0) el = std::fmin(std::fmin(2 * c, wb - c), pi);
        shch(c / 2, hs, hc);
    }
    rb = el; el = el / 4;
}

EXC_FN void rect_ratio(int kind1, double c, double hs, double hc, int kr, double y0, double yend,
                       double a, double b, double t, double &dd, double &de) {
    const int np = 10;
    double ct[np], st[np], f[np], f1[np];
    double el = pow2m(kr), x, sh, ch, sn, cn, rs, yv, y1, ff, ff1, q, g, g1, ser, ser1, pk, y00, ab, mt, sb, cb;
    const double *uu = exc_u8, *wq = exc_w8;
    if (kr > kmax) {
        for (int j = 0; j < np; ++j) { ct[j] = 1; st[j] = el * uu[j] / 2; }
    } else if (kind1 != 0) {
        const double *tc = kind1 < 0 ? exc_ch8h : exc_c8h, *ts = kind1 < 0 ? exc_sh8h : exc_s8h;
        for (int j = 0; j < np; ++j) { ct[j] = tc[5 * kr + hidx(j)]; st[j] = hsgn(j) * ts[5 * kr + hidx(j)]; }
    }
    switch (kind1) {
    case 0:
        for (int j = 0; j < np; ++j) {
            x = c + el * uu[j];
            ser = exc_cser[0]; ser1 = 0;
            for (int n = 1; n < 9; ++n) { ser1 = ser1 * x + ser; ser = ser * x + exc_cser[n]; }
            yv = y0 + b * x * ser; y1 = b * (ser + x * ser1);
            ff = exc_pser2[0]; ff1 = 0;
            for (int n = 1; n < np_ser; ++n) { ff1 = ff1 * x + ff; ff = ff * x + exc_pser2[n]; }
            q = std::sqrt(yv); g = yv * ff + a;
            f[j] = 1 / (q * g - t);
            f1[j] = -(.5 * y1 / q * g + q * (y1 * ff + yv * ff1)) * (f[j] * f[j]);
        }
        break;
    case 1: case 2: {
        bool us = (kind1 == 1) == (b >= 0);
        y00 = b >= 0 ? y0 : yend; ab = 2 * std::fabs(b);
        pk = kind1 == 2 ? -1 : 1;
        mt = kind1 == 1 ? 0.0 : pi;                   // m*pi or (m+1)*pi with m = 0
        for (int j = 0; j < np; ++j) {
            x = c + el * uu[j];
            sh = hs * ct[j] + hc * st[j]; ch = hc * ct[j] - hs * st[j];
            sn = 2 * sh * ch; cn = (ch - sh) * (ch + sh); rs = 1 / sn;
            double us2 = us ? sh : ch;
            yv = y00 + ab * (us2 * us2);
            y1 = pk * b * sn;
            ff = (pk * (x - sn * cn) + mt) * rsqrt2 * (rs * rs * rs);
            ff1 = (pk * sqrt2 - 3 * cn * ff) * rs;
            q = std::sqrt(yv); g = yv * ff + a;
            f[j] = 1 / (q * g - t);
            f1[j] = -(.5 * y1 / q * g + q * (y1 * ff + yv * ff1)) * (f[j] * f[j]);
        }
        break; }
    case -2:
        shch(yend, sb, cb);                           // yend carries wb
        for (int j = 0; j < np; ++j) {
            x = c + el * uu[j];
            sh = hs * ct[j] + hc * st[j]; ch = hc * ct[j] + hs * st[j];
            sn = sb * (1 + 2 * sh * sh) - cb * 2 * sh * ch; cn = cb * (1 + 2 * sh * sh) - sb * 2 * sh * ch;
            yv = 2 * b * (sb * ch - cb * sh) * sh; y1 = b * sn;
            fhyp_d(yend - x, sn, cn, ff, ff1);
            q = std::sqrt(yv); g = yv * ff + a;
            f[j] = 1 / (q * g - t);
            f1[j] = -(.5 * y1 / q * g + q * (y1 * ff - yv * ff1)) * (f[j] * f[j]);
        }
        break;
    default: {
        bool lw = b < 0;
        for (int j = 0; j < np; ++j) {
            x = c + el * uu[j];
            sh = hs * ct[j] + hc * st[j]; ch = hc * ct[j] + hs * st[j];
            sn = 2 * sh * ch; cn = 1 + 2 * (sh * sh); rs = 1 / sn;
            yv = y0 - 2 * b * (sh * sh); y1 = -b * sn;
            q = std::sqrt(yv);
            double rs2 = rs * rs, rs3 = rs2 * rs, rs4 = rs2 * rs2;
            if (lw) {
                g = ((y0 + b) * cn - b) * rs2 - yv * x * rs3;
                g1 = (y0 + b) * rs - 2 * ((y0 + b) * cn - b) * cn * rs3 - (y1 * x + yv) * rs3 + 3 * yv * x * cn * rs4;
                f[j] = 1 / (q * g * rsqrt2 - t);
                f1[j] = -(.5 * y1 / q * g + q * g1) * rsqrt2 * (f[j] * f[j]);
            } else {
                ff = (sn * cn - x) * rsqrt2 * rs3; ff1 = (sqrt2 - 3 * cn * ff) * rs;
                g = yv * ff + a;
                f[j] = 1 / (q * g - t);
                f1[j] = -(.5 * y1 / q * g + q * (y1 * ff + yv * ff1)) * (f[j] * f[j]);
            }
        }
        break; }
    }
    double fb = f[np - 2], fa = f[np - 1], q0 = 0, q1 = 0, q2 = 0, r0 = 0, r1 = 0, r2 = 0;
    for (int j = 0; j < np - 3; j += 2) {
        x = el * uu[j];
        q0 = q0 + wq[j] * f1[j]; q1 = q1 + wq[j] * (x * f1[j] + f[j]); q2 = q2 + wq[j] * x * (x * f1[j] + 2 * f[j]);
        x = el * uu[j + 1];
        r0 = r0 + wq[j + 1] * f1[j + 1]; r1 = r1 + wq[j + 1] * (x * f1[j + 1] + f[j + 1]);
        r2 = r2 + wq[j + 1] * x * (x * f1[j + 1] + 2 * f[j + 1]);
    }
    q0 = fb - fa - el * (q0 + r0);
    q1 = el * (fb + fa) - el * (q1 + r1);
    q2 = el * el * (fb - fa) - el * (q2 + r2);
    dd = q1 / q0;
    de = std::fabs(q2 / q0 - dd * dd);
    if (!(std::fabs(dd) < huge_)) {                   // exact node hit: that node is the root
        dd = huge_;
        for (int j = 0; j < np; ++j)
            if (!(std::fabs(f[j]) < huge_)) { dd = el * uu[j]; de = 0; }
    }
    if (!(de == de)) de = huge_;
}

EXC_FN void fast_accept(int kind1, double c, double dd, double &hs, double &hc, int &k, double &x, bool &ok) {
    ok = std::fabs(dd) < huge_;
    if (!ok) return;
    double cn = c + dd;
    if (kind1 == 0) {
        k = cn < 0 ? -1 : 1;
        x = std::sqrt(std::fabs(cn));
        if (k > 0) sincos_half(x, hs, hc); else shch(x / 2, hs, hc);
    } else {
        rotate(kind1, cn, dd / 2, hs, hc);
        k = kind1; x = cn;
    }
    ok = x > 0;
}

EXC_FN void fg_velocities(const double r1[3], const double r2[3], double yn, double scale, int d,
                          double v1[3], double v2[3]) {
    double y = yn * scale, n1 = std::sqrt(dot3(r1, r1)), n2 = std::sqrt(dot3(r2, r2));
    double ig = 1 / (d * std::sqrt(n1 * n2 + dot3(r1, r2)) * std::sqrt(y));
    double f1 = y / n1, f2 = y / n2;
    for (int k = 0; k < 3; ++k) {
        double s = r2[k] - r1[k];
        v1[k] = (s + f1 * r1[k]) * ig; v2[k] = (s - f2 * r2[k]) * ig;
    }
}

EXC_FN void velocities(const Frame &g, const double r1[3], const double r2[3], int direction, int k, double x,
                       double sw, double dw, bool fine, double v1[3], double v2[3]) {
    double p = sgn(1.0, g.a), c, s, y, dc, sh = sw, ch = dw, vr1, vr2, vt1, vt2, w;
    bool fg = fine;
    if (k > 0) { c = (ch - sh) * (ch + sh); w = 2 * sh * ch; s = w * w; }
    else { c = 1 + 2 * sh * sh; w = 2 * sh * ch; s = -(w * w); }
    if (k == 2) {
        c = -c;
        y = g.y0 + 2 * g.b * ch * ch;
        if (g.b < 0) y = g.yend - 2 * g.b * sh * sh;
    } else if (k > 0) {
        if (g.b >= 0) y = g.y0 + 2 * g.b * sh * sh; else y = g.yend - 2 * g.b * ch * ch;
    } else if (k == -2) {
        p = std::sqrt(g.y0 / (2 * g.b)); c = std::sqrt(1 + p * p);
        y = 2 * g.b * (2 * p * c * ch - (1 + 2 * p * p) * sh) * sh;
        dc = p * ch - c * sh; ch = c * ch - p * sh; sh = dc;
        c = 1 + 2 * sh * sh; w = 2 * sh * ch; s = -(w * w); fg = true;
        p = sgn(1.0, g.a);
    } else {
        y = g.y0 - 2 * g.b * sh * sh;
        if (g.b > 0 && y < g.y0 / 4) {
            fhyp_d(x, 2 * sh * ch, c, dc, s); dc = y * dc;
            if (dc < g.a) {
                y = g.t / (dc + g.a); y = y * y;
                sh = std::sqrt((g.y0 - y) / (2 * g.b)); ch = std::sqrt(1 + sh * sh); c = 1 + 2 * sh * sh;
                w = 2 * sh * ch; s = -(w * w);
                fg = true;
            }
        }
    }
    dc = p * g.cs - c;
    if (p * g.cs * c > 0 && std::fabs(p * g.cs + c) > 0.5) dc = (s - g.ss * g.ss) / (p * g.cs + c);
    if (fg && g.cs > .3) {
        fg_velocities(r1, r2, y, g.scale, direction, v1, v2);
    } else {
        y = sqrt2 / std::sqrt(y * g.scale);
        dc = y * dc; p = y * p * g.cs * g.q; s = y * g.ss;
        sh = 1 / g.s1; ch = 1 / g.s2;
        vr1 = dc - p * sh; vt1 = s * g.s2 * sh;
        vr2 = -(dc + p * ch); vt2 = s * g.s1 * ch;
        for (int i = 0; i < 3; ++i) {
            v1[i] = vr1 * g.e1[i] + vt1 * g.t1v[i];
            v2[i] = vr2 * g.e2[i] + vt2 * g.t2v[i];
        }
    }
}

// excali_one for nrev = 0 (normalized inputs, mu = 1): status 0, 2 (low
// confidence); negative on invalid input or failure; -6 for nrev /= 0 (not ported).
EXC_FN void excali_zero_rev(const double r1[3], const double r2[3], double tof, int direction,
                            double v1[3], double v2[3], int &status, int &half, int nrev = 0,
                            const double *plane_normal = nullptr) {
    const double nan = std::nan("");
    Frame g;
    double x, c = 0, hs, hc, wb, dd, de, rb, el;
    int k, kd = 0, kr;
    bool fast, nearp, fin;
    for (int i = 0; i < 3; ++i) { v1[i] = nan; v2[i] = nan; }
    half = 0; status = -1;
    if ((direction != 1 && direction != -1) || nrev < 0 || !(tof > 0) || !std::isfinite(tof)) return;
    if (nrev != 0) { status = -6; return; }
    geometry(r1, r2, tof, direction, g, status, half, plane_normal);
    if (status != -4) return;
    rb = 0; kr = 0; hs = 0; hc = 1; dd = 0; de = 0; wb = 0;
    zero_rev_start(g, k, x, wb, fin, nearp);
    fast = x > 0 || k == 0;
    if (fast && k != 0) {
        fast_coordinate(k, x, g.y0, g.b, wb, kd, c, el, rb, hs, hc);
        fast = el > 0; kr = kexp(el);
        if (fast && !fin) rect_ratio(kd, c, hs, hc, kr, g.y0, kd == -2 ? wb : g.yend, g.a, g.b, g.t, dd, de);
        if (fast) fast_accept(kd, c, dd, hs, hc, k, x, fast);
    }
    if (fast) {
        status = de > 1e-8 * rb * pow2m(kr) ? 2 : 0;
        if (k == 0) { k = 1; x = 0; hs = 0; hc = 1; }
    }
    if (status >= 0 && status <= 2) {
        velocities(g, r1, r2, direction, k, x, hs, hc, nearp, v1, v2);
        if (!finite3(v1) || !finite3(v2)) status = -4;
    }
    if (status < 0 || status > 2) for (int i = 0; i < 3; ++i) { v1[i] = nan; v2[i] = nan; }
}

}  // namespace exc
