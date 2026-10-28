"""Generate EXCALI's compiled-in starter tables (src/excali_tables.inc).

usage: python3 tools/gen_tables.py [N] [OUT]   (shipped: N=32)

A zero-revolution root depends only on the Lancaster-Blanchard parameter lam
and on the time relative to the parabolic time tp. All grids are uniform in
tau = lam/(1+sqrt(1-lam**2)).

Ellipse (t > tp), root v (half the eccentric-anomaly change), s = (tp/t)**(1/3),
|tau| <= 0.9: tab(:,:,1) = (pi - v)/s for s in [0, 1/2] and tab(:,:,2) =
v**2/(1 - s) for s in [1/2, 1], smooth through t -> infinity and t -> tp.

Hyperbola (t < tp), root w (half the hyperbolic-anomaly change), u = t/tp,
on sqrt|tau| uniform in [0.02, sqrt(0.9)]: short way (u in [0.1, 1])
tabh(:,:,1) = w**2/((1-u)*wb**2) with wb the branch point of sqrt(y); long way
(u = p**3, p in [0.05, 1]) tabh(:,:,2) = w**2/((1-u)*(1-log(u))**2); very short
short-way hyperbolas (u in [0, 0.1]) tabh(:,:,3) = xi/(wb*u**2) with xi = wb - w
the distance from the branch point.

Multiple revolutions, slices in mu = 1/m (m in MR_M, the last ~infinity),
interpolated in mu at run time, |tau| <= TM:
  mr_c(:,k) = cbrt(t_min/m)/sqrt(s_L) and mr_s(:,k) = v* (the minimum-time
  location) on a fine tau grid of MR_NC points;
  mr_v(:,:,k) = branch root v on tau x sigma, sigma = branch*sqrt(1-(t_min/t)**(1/3)).
The starters only need to be within a few percent of the singularity
distance, so their grids are coarse. mr_c is not a starter: it decides
feasibility outside a 1e-6 margin, so it needs MR_NC >= 128 and the m = 5
slice (worst relative error 1.04e-7 for m = 1..40 and sampled m to 1000,
including m between slices; 4.4e-7 at m = 5 without that slice).

Truth: bisection on EXCALI's T for the canonical geometry r1 = r2 = 1, with
the parabolic series where the closed form cancels.
"""
import sys
import numpy as np

S2 = np.sqrt(2.0); PI = np.pi; TM = 0.9
PSER = np.array([115665628927/2601853536522240000, 3821869001/10751460894720000,
                 133541/48117888000, 763967/36324288000, 1393481/9081072000,
                 1181/1108800, 29/4200, 17/420, 1/5, 2/3])


def geom(lam):
    """y0, yend, a, b and tp of the canonical geometry r1 = r2 = 1 with this lam."""
    d = np.where(lam >= 0, 1.0, -1.0)
    ss = (1 - lam * lam) / (1 + lam * lam); cs = np.sqrt(np.maximum(0, 1 - ss * ss))
    a = d * S2 * cs; b = S2 * a
    yend = 2 * ss * ss / (1 + cs)
    y0 = np.where(d > 0, yend, 2 * (1 + cs))
    tp = np.sqrt(y0) * (y0 * S2 / 3 + a)
    return a, b, y0, yend, tp


def T_v(v, a, b, y0, yend, m=0):
    sh, ch = np.sin(v / 2), np.cos(v / 2)
    y = np.where(b < 0, yend - 2 * b * ch * ch, y0 + 2 * b * sh * sh)
    sv, cv = np.sin(v), np.cos(v)
    x = v * v; ser = np.zeros_like(v)
    for p in PSER:
        ser = ser * x + p
    F = np.where((v < .3) & (m == 0), ser / S2, (v + m * PI - sv * cv) / (S2 * sv ** 3))
    return np.sqrt(y) * (y * F + a)


def root_v(lam, s):
    a, b, y0, yend, tp = geom(lam); t = tp / s ** 3
    lo = np.zeros_like(lam); hi = np.full(lam.shape, PI)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f = T_v(mid, a, b, y0, yend) < t
        lo = np.where(f, mid, lo); hi = np.where(f, hi, mid)
    return 0.5 * (lo + hi)


def T_w(w, a, b, y0):
    """Hyperbolic T(w); cancellation-free long-way form for b < 0, w >= 0.3."""
    sh = np.sinh(w / 2); s = np.sinh(w); c = np.cosh(w)
    y = y0 - 2 * b * sh * sh
    q = np.sqrt(np.maximum(y, 0))
    X = -w * w; ser = np.zeros_like(w)
    for p in PSER:
        ser = ser * X + p
    sb = np.where(w < .3, 1, s)
    F = np.where(w < .3, ser / S2, (s * c - w) / (S2 * sb ** 3))
    Tlw = q * (((y0 + b) * c - b) / sb ** 2 - y * w / sb ** 3) / S2
    return np.where((b < 0) & (w >= .3), Tlw, q * (y * F + a))


def branch_w(b, y0):
    return 2 * np.arcsinh(np.sqrt(y0 / (2 * np.abs(b))))


def root_w(lam, u):
    a, b, y0, yend, tp = geom(lam); t = u * tp
    if np.all(b > 0) and np.all(u <= PX):          # in xi = wb - w, exact near the branch point
        wb = branch_w(b, y0); lo = np.zeros_like(lam); hi = wb.copy()
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            w = wb - mid; sh = np.sinh(mid / 2); y = 2 * b * np.sinh(wb - mid / 2) * sh
            s = np.sinh(w); c = np.cosh(w); X = -w * w; ser = np.zeros_like(w)
            for p in PSER:
                ser = ser * X + p
            F = np.where(w < .3, ser / S2, (s * c - w) / (S2 * np.where(w < .3, 1, s) ** 3))
            f = np.sqrt(np.maximum(y, 0)) * (y * F + a) < t     # T increases with xi
            lo = np.where(f, mid, lo); hi = np.where(f, hi, mid)
        return wb - 0.5 * (lo + hi)
    lo = np.zeros_like(lam); hi = np.where(b > 0, branch_w(b, y0), 80.0)
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        f = T_w(mid, a, b, y0) > t            # T decreases in w
        lo = np.where(f, mid, lo); hi = np.where(f, hi, mid)
    return 0.5 * (lo + hi)


SH = (0.02, 0.9 ** 0.5); PS = (0.1, 1.0); PL = (0.05, 1.0); PX = 0.1


def hyperbolic(n):
    """Hyperbolic slices on sqrt|tau| uniform in SH (fine toward 180 degrees)."""
    s = np.linspace(SH[0], SH[1], n)
    out = []
    for side, (p0, p1) in ((1, PS), (-1, PL)):
        tau = side * s * s; lam = 2 * tau / (1 + tau * tau)
        L, P = np.meshgrid(lam, np.linspace(p0, p1, n), indexing='ij')
        u = np.minimum(P if side > 0 else P ** 3, 1 - 1e-6)
        w = root_w(L, u)
        a, b, y0, yend, tp = geom(L)
        if side > 0:
            out.append(w * w / ((1 - u) * branch_w(b, y0) ** 2))
        else:
            out.append(w * w / ((1 - u) * (1 - np.log(u)) ** 2))
    # very short short-way hyperbolas: xi = wb - w on u = PX*p
    tau = s * s; lam = 2 * tau / (1 + tau * tau)
    L, P = np.meshgrid(lam, np.linspace(0, 1, n), indexing='ij')
    u = np.maximum(PX * P, 1e-7)
    a, b, y0, yend, tp = geom(L); wb = branch_w(b, y0)
    out.append((wb - root_w(L, u)) / (wb * u * u))
    return out


MR_M = [1, 2, 3, 4, 5, 6, 8, 12, 16, 32, 1000000]; MR_NT = 16; MR_NS = 16; MR_NC = 128


def tmin_v(lam, m):
    """Minimum time and its location v* (bisection on dT/dv; T is unimodal)."""
    a, b, y0, yend, tp = geom(lam)
    lo = np.full(np.shape(lam), 1e-3); hi = np.full(np.shape(lam), PI - 1e-3)
    h = 1e-6
    for _ in range(60):
        mid = .5 * (lo + hi)
        f = (T_v(mid + h, a, b, y0, yend, m) - T_v(mid - h, a, b, y0, yend, m)) < 0
        lo = np.where(f, mid, lo); hi = np.where(f, hi, mid)
    vs = .5 * (lo + hi)
    return T_v(vs, a, b, y0, yend, m), vs


def branch_v(lam, sig, m):
    a, b, y0, yend, tp = geom(lam); tm, vs = tmin_v(lam, m)
    t = tm / (1 - sig * sig) ** 3
    lower = sig < 0
    lo = np.where(lower, 1e-12, vs); hi = np.where(lower, vs, PI - 1e-12)
    for _ in range(80):
        mid = .5 * (lo + hi); Tv = T_v(mid, a, b, y0, yend, m)
        right = np.where(lower, Tv > t, Tv < t)
        lo = np.where(right, mid, lo); hi = np.where(right, hi, mid)
    return .5 * (lo + hi)


def multi_rev(n):
    tg = np.linspace(-TM, TM, n); lg = 2 * tg / (1 + tg * tg); sg = np.linspace(-1, 1, MR_NS)
    L, S = np.meshgrid(lg, np.clip(sg, -1 + 1e-9, 1 - 1e-9), indexing='ij')
    tc = np.linspace(-TM, TM, MR_NC); lc = 2 * tc / (1 + tc * tc)
    V = np.zeros((n, MR_NS, len(MR_M))); C = np.zeros((MR_NC, len(MR_M))); VS = np.zeros_like(C)
    for k, m in enumerate(MR_M):
        V[:, :, k] = branch_v(L, S, m)
        tm, vs = tmin_v(lc, m)
        # T scales as s_L**1.5 (canonical s_L = 2/(1+lam**2)).
        C[:, k] = np.cbrt(tm / m) / np.sqrt(2 / (1 + lc * lc)); VS[:, k] = vs
    return V, C, VS


def zero_rev(n):
    tau = np.linspace(-TM, TM, n); lam = 2 * tau / (1 + tau * tau)
    su = np.linspace(0, .5, n); sl = np.linspace(.5, 1, n)
    L, SU = np.meshgrid(lam, su, indexing='ij'); _, SL = np.meshgrid(lam, sl, indexing='ij')
    SU = np.maximum(SU, 1e-7); SL = np.minimum(SL, 1 - 1e-7)
    A = (PI - root_v(L, SU)) / SU
    B = root_v(L, SL) ** 2 / (1 - SL)
    return A, B


def lagrange4(p):
    return np.stack([-(p - 1) * (p - 2) * (p - 3) / 6, p * (p - 2) * (p - 3) / 2,
                     -p * (p - 1) * (p - 3) / 2, p * (p - 1) * (p - 2) / 6])


def interp(M, n, u_tau, u_s):
    i = np.clip(np.floor(u_tau).astype(int) - 1, 0, n - 4); wt = lagrange4(u_tau - i)
    j = np.clip(np.floor(u_s).astype(int) - 1, 0, n - 4); ws = lagrange4(u_s - j)
    acc = np.zeros_like(u_tau)
    for a in range(4):
        for b in range(4):
            acc += wt[a] * ws[b] * M[i + a, j + b]
    return acc


def estimate(A, B, n, lam, s):
    tau = lam / (1 + np.sqrt(1 - lam * lam))
    ut = (tau + TM) * ((n - 1) / (2 * TM))
    vA = PI - s * interp(A, n, ut, s * 2 * (n - 1))
    vB = np.sqrt(np.maximum(0, (1 - s) * interp(B, n, ut, (s - .5) * 2 * (n - 1))))
    return np.where(s <= .5, vA, vB)


def arr(f, name, M):
    """Declare a protected array and initialize it with DATA statements of
    whole columns (at most 1000 values, under the standard's 255
    continuation lines)."""
    shp = M.shape
    f.write(f"   real(dp),protected::{name}({','.join(str(x) for x in shp)})\n")
    n1 = shp[0]; cols = [()] if len(shp) == 1 else list(np.ndindex(*shp[1:][::-1]))
    cols = [c[::-1] for c in cols]                      # Fortran column order
    per = max(1, 1000 // n1)
    groups = []
    for c in cols:
        if groups and len(groups[-1]) < per and groups[-1][-1][1:] == c[1:]:
            groups[-1].append(c)
        else:
            groups.append([c])
    for grp in groups:
        if len(shp) == 1:
            target = f"{name}(1:{n1})"
        else:
            # consecutive first-trailing indices with the rest fixed
            lead = [c[0] for c in grp]; rest = grp[0][1:]
            assert all(c[1:] == rest for c in grp) and lead == list(range(lead[0], lead[0] + len(lead)))
            idx = ",".join(str(r + 1) for r in rest)
            target = f"(({name}(ti_,tj_{',' + idx if idx else ''}),ti_=1,{n1}),tj_={lead[0] + 1},{lead[-1] + 1})"
        vals = np.concatenate([M[(slice(None),) + c] for c in grp]) if len(shp) > 1 else M
        f.write(f"   data {target}/ &\n")
        f.write(", &\n".join("      " + ",".join(f"{x:.17e}_dp" for x in vals[k:k + 4]) for k in range(0, len(vals), 4)))
        f.write("/\n")


def write(out, n, A, B, HS, HL, HX, V, C, VS):
    def block(f, name, mats):
        arr(f, name, np.stack(mats, axis=-1))
    with open(out, 'w') as f:
        f.write(f"   ! Generated by tools/gen_tables.py {n}; do not edit.\n")
        f.write(f"   ! Zero-revolution root tables on tau = lam/(1+sqrt(1-lam**2)).\n")
        f.write(f"   ! Ellipse, |tau| <= {TM}, s = (tp/t)**(1/3): tab(:,:,1) = (pi-v)/s on [0,1/2],\n")
        f.write("   ! tab(:,:,2) = v**2/(1-s) on [1/2,1].\n")
        f.write(f"   ! Hyperbola, sqrt|tau| in [{SH[0]},{SH[1]:.6f}], u = t/tp: tabh(:,:,1) = w**2/((1-u)*wb**2)\n")
        f.write(f"   ! (short way, u in [{PS[0]},1]); tabh(:,:,2) = w**2/((1-u)*(1-log(u))**2)\n")
        f.write(f"   ! (long way, u = p**3, p in [{PL[0]},1]); tabh(:,:,3) = xi/(wb*u**2), xi = wb-w\n")
        f.write(f"   ! (short way, u in [0,{PX}]).\n")
        f.write(f"   integer,parameter::tab_n={n}\n   real(dp),parameter::tab_tmax={TM}_dp\n")
        f.write("   integer,private::ti_,tj_\n")
        f.write(f"   real(dp),parameter::tabh_s0={SH[0]}_dp,tabh_s1={SH[1]:.17e}_dp,tabh_ps={PS[0]}_dp,tabh_pl={PL[0]}_dp,tabh_px={PX}_dp\n")
        block(f, "tab", (A, B))
        block(f, "tabh", (HS, HL, HX))
        f.write(f"   ! Multiple revolutions (slices in mu = 1/m, last ~infinity): mr_v(tau,sigma,k) branch root,\n")
        f.write(f"   ! mr_c(tau,k) = cbrt(t_min/m)/sqrt(s_L), mr_s(tau,k) = v*; |tau| <= {TM}.\n")
        f.write(f"   integer,parameter::mr_nt={V.shape[0]},mr_ns={V.shape[1]},mr_nm={len(MR_M)},mr_nc={MR_NC}\n")
        f.write("   real(dp),parameter::mr_mu(mr_nm)=[" + ",".join(f"{1 / m:.17e}_dp" for m in MR_M) + "]\n")
        f.write("   integer,parameter::mr_m(mr_nm)=[" + ",".join(str(m) for m in MR_M) + "]\n")
        arr(f, "mr_v", V); arr(f, "mr_c", C); arr(f, "mr_s", VS)

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    out = sys.argv[2] if len(sys.argv) > 2 else "src/excali_tables.inc"
    A, B = zero_rev(n)
    HS, HL, HX = hyperbolic(n)
    V, C, VS = multi_rev(MR_NT)
    write(out, n, A, B, HS, HL, HX, V, C, VS)
