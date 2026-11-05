"""Generate EXCALI's compiled-in starter tables (src/excali_tables.inc).

usage: python3 tools/gen_tables.py [N] [OUT]   (shipped: N=32)

A zero-revolution root depends only on the Lancaster-Blanchard parameter lam
and on the time relative to the parabolic time tp. All grids are uniform in
tau = lam/(1+sqrt(1-lam**2)).

Short-way ellipses toward tau -> 1 (tau in PP): tabp, split at v = pi/2; see
short_way().

Ellipse (t > tp), root v (half the eccentric-anomaly change), s = (tp/t)**(1/3),
|tau| <= 0.9 (the solver uses it below tau = 0.75): tab(:,:,1) = (pi - v)/(s*ck) for s in [0, 1/2], ck = (K/tp)**(1/3)
with the pole asymptote T ~ K/(pi-v)**3, K = pi*y(pi)**1.5/sqrt2, and
tab(:,:,2) = v**2*T0'/(t - tp) for s in [1/2, 1], T0' = dT/d(v**2) at the
parabola. Both tend to 1 at their ends (t -> infinity, t -> tp), and the
normalization absorbs the sqrt(y) branch point's approach to the parabola as
|tau| grows, so the table is about three times more accurate there.

Hyperbola (t < tp), root w (half the hyperbolic-anomaly change), u = t/tp in
[UH, 1], on sqrt|tau| uniform in [0, sqrt(TOP)]: short way
tabh(:,:,1) = w**2/((1-u)*wn**2), wn**2 = L**2*(1-exp(-wb**2/L**2)), L = 1-log(u),
with wb the branch point of sqrt(y): wn -> wb as tau -> 1, and wn -> L as
tau -> 0 (wb ~ log(1/tau)), where the two ways meet, with corrections that
vanish faster than any power of tau (a rational blend leaves a 1/log(tau)
term the grid cannot follow: 1.3% instead of 0.06%); long way tabh(:,:,2) = w**2/((1-u)*L**2). Below UH
the solver uses the large-w closed form instead.

Multiple revolutions, slices in mu = 1/m (m in MR_M, the last ~infinity),
interpolated in mu at run time, |tau| <= TM:
  mr_c(:,k) = cbrt(t_min/m)/sqrt(s_L) and mr_s(:,k) = v* (the minimum-time
  location) on a fine tau grid of MR_NC points;
  mr_q(:,:,h,k) = branch root normalized by its two anchors on tau x |sigma|,
  sigma = branch*sqrt(1-(t_min/t)**(1/3)), one half per branch (h=1: sigma<0,
  h=2: sigma>0): q = (v - v*)/(sigma*(P - v*)), P = 0 or pi the pole on the
  branch's side. q is smooth on each half, tends to sqrt(6 t_min/T''(v*))/|P-v*|
  at the fold and to 1 at the pole, where v = v* + sigma*(P-v*)*q hits P exactly
  whatever the error of the interpolated v*. With 32 tau points and 16 |sigma|
  points it is within 2.5% of |v - v*| (or the pole distance) for |sigma| >=
  1e-3, so the rectangle needs no certificate (the old 16x16 table of v itself
  was off by up to 8 times |v - v*| near the fold). mr_c decides
feasibility outside a 1e-6 margin, so it needs MR_NC >= 128 and the m = 5
slice (worst relative error 1.04e-7 for m = 1..40 and sampled m to 1000,
including m between slices; 4.4e-7 at m = 5 when that slice is left out).

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
    lo = np.zeros_like(lam); hi = np.where(b > 0, branch_w(b, y0), 80.0)
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        f = T_w(mid, a, b, y0) > t            # T decreases in w
        lo = np.where(f, mid, lo); hi = np.where(f, hi, mid)
    return 0.5 * (lo + hi)


UH = 0.1; TOP = 1 - 1e-6           # hyperbolic table: u = t/tp in [UH, 1], |tau| <= TOP


def hyperbolic(n):
    """Hyperbolic slices on sqrt|tau| uniform in [0, sqrt(TOP)], u uniform in [UH, 1]."""
    s = np.linspace(0, np.sqrt(TOP), n)
    out = []
    for side in (1, -1):
        tau = side * s * s; lam = 2 * tau / (1 + tau * tau)
        L, P = np.meshgrid(lam, np.linspace(UH, 1, n), indexing='ij')
        u = np.minimum(P, 1 - 1e-6)
        w = root_w(L, u)
        a, b, y0, yend, tp = geom(L)
        lg = (1 - np.log(u)) ** 2
        if side > 0:
            with np.errstate(divide='ignore'):
                wb2 = branch_w(b, y0) ** 2               # infinite at tau = 0: wn = L
            out.append(w * w / ((1 - u) * lg * -np.expm1(-wb2 / lg)))
        else:
            out.append(w * w / ((1 - u) * lg))
    return out


MR_M = [1, 2, 3, 4, 5, 6, 8, 12, 16, 32, 1000000]; MR_NT = 32; MR_NS = 16; MR_NC = 128


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
    tg = np.linspace(-TM, TM, n); lg = 2 * tg / (1 + tg * tg); sg = np.linspace(0, 1, MR_NS)
    tc = np.linspace(-TM, TM, MR_NC); lc = 2 * tc / (1 + tc * tc)
    V = np.zeros((n, MR_NS, 2, len(MR_M))); C = np.zeros((MR_NC, len(MR_M))); VS = np.zeros_like(C)
    for k, m in enumerate(MR_M):
        _, vg = tmin_v(lg, m)
        for h, br in enumerate((-1, 1)):
            S = br * np.clip(sg, 1e-7, 1 - 1e-9)
            L, SS = np.meshgrid(lg, S, indexing='ij')
            P = np.pi if br > 0 else 0.0
            V[:, :, h, k] = (branch_v(L, SS, m) - vg[:, None]) / (SS * (P - vg[:, None]))
        tm, vs = tmin_v(lc, m)
        # T scales as s_L**1.5 (canonical s_L = 2/(1+lam**2)).
        C[:, k] = np.cbrt(tm / m) / np.sqrt(2 / (1 + lc * lc)); VS[:, k] = vs
    return V, C, VS


def zero_rev(n):
    tau = np.linspace(-TM, TM, n); lam = 2 * tau / (1 + tau * tau)
    su = np.linspace(0, .5, n); sl = np.linspace(.5, 1, n)
    L, SU = np.meshgrid(lam, su, indexing='ij'); _, SL = np.meshgrid(lam, sl, indexing='ij')
    SU = np.maximum(SU, 1e-7); SL = np.minimum(SL, 1 - 1e-7)
    ck, dt = ell_scales(L)
    A = (PI - root_v(L, SU)) / (SU * ck)
    B = root_v(L, SL) ** 2 / ((SL ** -3 - 1) * dt)
    return A, B


PP = (0.75, 1 - 1e-6)                # tau range of the short-way table tabp


def zp(tau):
    """Grid coordinate of tabp: 1 - sqrt(1 - tau), fine as tau -> 1."""
    return 1 - np.sqrt(1 - tau)


def half_scales(lam):
    """tp, T_half = T(pi/2) and K = pi*y(pi)**1.5/sqrt2 for the canonical geometry."""
    a, b, y0, yend, tp = geom(lam)
    ypi = np.where(b >= 0, y0 + 2 * b, yend)
    yh = np.where(b >= 0, y0 + b, yend - b)
    return tp, np.sqrt(yh) * (yh * PI / (2 * S2) + a), PI * ypi ** 1.5 / S2


def short_way(n):
    """Short-way ellipses toward tau -> 1, split at v = pi/2 (T_half):
    tabp(:,:,1) = (pi-v)/(s1*(K/T_half)**(1/3)) on s1 = (T_half/t)**(1/3), and
    tabp(:,:,2) = v**2/(s2*(s2+2*e0)) on s2 = (t-tp)/(T_half-tp), e0 = tp/(T_half-tp),
    which absorbs the sqrt(y) branch point's approach to the parabola."""
    tau = 1 - (1 - np.linspace(zp(PP[0]), zp(PP[1]), n)) ** 2; lam = 2 * tau / (1 + tau * tau)
    L, S = np.meshgrid(lam, np.linspace(0, 1, n), indexing='ij')
    tp, th, K = half_scales(L); a, b, y0, yend, _ = geom(L)
    s1 = np.maximum(S, 1e-7)
    A = (PI - root_t(L, th / s1 ** 3)) / (s1 * np.cbrt(K / th))
    s2 = np.clip(S, 1e-7, 1 - 1e-12); e0 = tp / (th - tp)
    B = root_t(L, tp + s2 * (th - tp)) ** 2 / (s2 * (s2 + 2 * e0))
    return A, B


NB = 64


def P_body(d):
    """T/|b|**1.5 near the long-way limit y_end -> 0, d = pi - v."""
    return (PI - d + np.sin(d) * np.cos(d)) / (4 * np.cos(d / 2) ** 3) - np.sin(d / 2)


def long_body(nb=NB):
    """Long-way limit tau -> -1: d = pi - v from P_body(d) = t/|b|**1.5 in
    u = (t/|b|**1.5 - 1/3)/(pi/4 - 1/3): tabb(:,1) = v**2/u on u in [0,1/2],
    tabb(:,2) = d/(1-u) on [1/2,1] (P_body is 1/3 at the parabola, pi/4 at the pole)."""
    ul = np.clip(np.linspace(0, .5, nb), 1e-6, None); uh = np.clip(np.linspace(.5, 1, nb), None, 1 - 1e-6)
    def inv(th):
        lo = np.zeros_like(th); hi = np.full(th.shape, PI - 1e-12)
        for _ in range(90):
            m = .5 * (lo + hi); f = P_body(m) > th
            lo = np.where(f, m, lo); hi = np.where(f, hi, m)
        return .5 * (lo + hi)
    vl = PI - inv(1 / 3 + ul * (PI / 4 - 1 / 3)); dh = inv(1 / 3 + uh * (PI / 4 - 1 / 3))
    return np.stack([vl ** 2 / ul, dh / (1 - uh)], axis=-1)


def root_t(lam, t):
    """Elliptic root v of T(v) = t by bisection."""
    a, b, y0, yend, tp = geom(lam)
    lo = np.zeros_like(lam); hi = np.full(lam.shape, PI)
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        f = T_v(mid, a, b, y0, yend) < t
        lo = np.where(f, mid, lo); hi = np.where(f, hi, mid)
    return 0.5 * (lo + hi)


def ell_scales(lam):
    """ck = (K/tp)**(1/3) (pole asymptote) and tp/T0' (parabola), both scale-free."""
    a, b, y0, yend, tp = geom(lam)
    ypi = np.where(b >= 0, y0 + 2 * b, yend)
    tv2 = y0 * np.sqrt(y0) / (5 * S2) + b * (np.sqrt(y0) / S2 + a / (2 * np.sqrt(y0))) / 2
    return np.cbrt(PI * ypi ** 1.5 / (S2 * tp)), tp / tv2


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
    ck, dt = ell_scales(lam)                    # the solver computes the one it needs
    vA = PI - s * ck * interp(A, n, ut, s * 2 * (n - 1))
    vB = np.sqrt(np.maximum(0, (s ** -3 - 1) * dt * interp(B, n, ut, (s - .5) * 2 * (n - 1))))
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


def write(out, n, A, B, PA, PB, LB, HS, HL, V, C, VS):
    def block(f, name, mats):
        arr(f, name, np.stack(mats, axis=-1))
    with open(out, 'w') as f:
        f.write(f"   ! Generated by tools/gen_tables.py {n}; do not edit.\n")
        f.write(f"   ! Zero-revolution root tables on tau = lam/(1+sqrt(1-lam**2)).\n")
        f.write(f"   ! Ellipse, |tau| <= {TM}, s = (tp/t)**(1/3): tab(:,:,1) = (pi-v)/(s*ck) on [0,1/2],\n")
        f.write("   ! ck = (pi*y(pi)**1.5/(sqrt2*tp))**(1/3); tab(:,:,2) = v**2*T0'/(t-tp) on [1/2,1],\n")
        f.write("   ! T0' = dT/d(v**2) at the parabola.\n")
        f.write(f"   ! Hyperbola, sqrt|tau| in [0,sqrt({TOP})], u = t/tp in [{UH},1]: tabh(:,:,1) = w**2/((1-u)*wn**2),\n")
        f.write("   ! wn**2 = L**2*(1-exp(-wb**2/L**2)), L = 1-log(u) (short way); tabh(:,:,2) = w**2/((1-u)*L**2) (long way).\n")
        f.write(f"   integer,parameter::tab_n={n}\n   real(dp),parameter::tab_tmax={TM}_dp\n")
        f.write("   integer,private::ti_,tj_\n")
        f.write(f"   real(dp),parameter::tabh_s1={np.sqrt(TOP):.17e}_dp,tabh_u0={UH}_dp\n")
        block(f, "tab", (A, B))
        f.write(f"   ! Short-way ellipses, tau in [{PP[0]},1-{1-PP[1]:.0e}], grid uniform in 1-sqrt(1-tau), split\n")
        f.write("   ! at v=pi/2: tabp(:,:,1)=(pi-v)/(s1*(K/Th)**(1/3)), s1=(Th/t)**(1/3), Th=T(pi/2);\n")
        f.write("   ! tabp(:,:,2)=v**2/(s2*(s2+2*e0)), s2=(t-tp)/(Th-tp), e0=tp/(Th-tp).\n")
        f.write(f"   real(dp),parameter::tabp_z0={zp(PP[0]):.17e}_dp,tabp_z1={zp(PP[1]):.17e}_dp\n")
        block(f, "tabp", (PA, PB))
        f.write("   ! Long-way limit tau -> -1, body d=pi-v from T/|b|**1.5=th, u=(th-1/3)/(pi/4-1/3):\n")
        f.write("   ! tabb(:,1)=v**2/u on [0,1/2], tabb(:,2)=d/(1-u) on [1/2,1].\n")
        f.write(f"   integer,parameter::tabb_n={LB.shape[0]}\n")
        arr(f, "tabb", LB)
        block(f, "tabh", (HS, HL))
        f.write(f"   ! Multiple revolutions (slices in mu = 1/m, last ~infinity): mr_q(tau,|sigma|,branch,k) = (v-v*)/(sigma*(P-v*)),\n")
        f.write(f"   ! mr_c(tau,k) = cbrt(t_min/m)/sqrt(s_L), mr_s(tau,k) = v*; |tau| <= {TM}.\n")
        f.write(f"   integer,parameter::mr_nt={V.shape[0]},mr_ns={V.shape[1]},mr_nm={len(MR_M)},mr_nc={MR_NC}\n")
        f.write("   real(dp),parameter::mr_mu(mr_nm)=[" + ",".join(f"{1 / m:.17e}_dp" for m in MR_M) + "]\n")
        f.write("   integer,parameter::mr_m(mr_nm)=[" + ",".join(str(m) for m in MR_M) + "]\n")
        arr(f, "mr_q", V); arr(f, "mr_c", C); arr(f, "mr_s", VS)

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    out = sys.argv[2] if len(sys.argv) > 2 else "src/excali_tables.inc"
    A, B = zero_rev(n)
    PA, PB = short_way(n)
    LB = long_body()
    HS, HL = hyperbolic(n)
    V, C, VS = multi_rev(MR_NT)
    write(out, n, A, B, PA, PB, LB, HS, HL, V, C, VS)
