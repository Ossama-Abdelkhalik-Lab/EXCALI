"""Generate EXCALI's compiled-in starter tables (src/excali_tables.inc).

usage: python3 tools/gen_tables.py [N] [OUT]   (shipped: N=32)

A zero-revolution root depends only on the Lancaster-Blanchard parameter lam
and on the time relative to the parabolic time tp. All grids are uniform in
tau = lam/(1+sqrt(1-lam**2)).

Ellipse (t > tp), root v (half the eccentric-anomaly change): one table tabe
for tau in [TE0, TE1] on a grid uniform in z = sign(tau)*(1 - sqrt(1 - |tau|))
(fine toward both ends), split at v = pi/2 (T_half = T(pi/2)); see ellipse().
Below TE0 (transfers near 360 degrees) the solver uses the closed-form
long-way starter with a polynomial body fit.

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
  h=2: sigma>0): q = (v - v*)/(|sigma|*(P - v*)), P = 0 or pi the pole on the
  branch's side. q is smooth on each half, tends to sqrt(6 t_min/T''(v*))/|P-v*|
  at the fold and to 1 at the pole, where v = v* + |sigma|*(P-v*)*q hits P exactly
  whatever the error of the interpolated v*. The table stores
  g = (q - 1)/(1 - |sigma|), finite at the pole, so the pole distance
  (1 - |sigma| q)|P - v*| = (1-|sigma|)(1 - |sigma| g)|P - v*| with
  1 - |sigma| = r/(1+|sigma|), r = (t_min/t)**(1/3), keeps full precision as
  t -> infinity (storing q itself leaves a 1e-9 offset at the pole). Its grid
  is 16 tau points clustered toward |tau| = TM (tau = TM*sin(1.45x)/sin(1.45),
  x uniform) by 12 points uniform in sqrt|sigma|: within 1.4% of |v - v*| (or
  the pole distance) for |sigma| >= 1e-2 in 4,000 random cases, where the
  uniform 32 x 16 grid it replaced reached 2.8% (the old 16x16 table of v
  itself was off by up to 8 times |v - v*| near the fold). mr_c decides
feasibility outside a 1e-6 margin, so it needs MR_NC >= 128 and the m = 5
slice (worst relative error 1.04e-7 for m = 1..40 and sampled m to 1000,
including m between slices; 4.4e-7 at m = 5 when that slice is left out).

Body limits of multiple revolutions off the tables (|tau| > 0.9): with
mu = 1/M (short way, M = m*pi) or -1/M (long way, M = (m+1)*pi), the large-|b|
time is T = M*|b|**1.5*Q(x)/4, Q(x) = (1 + mu*(x + sin x))/cos(x/2)**3, x = v
(short way) or pi - v (long way). Q increases to infinity at pi; on the long way
it has a minimum at x_b ~ 8|mu|/3, where db/(8|mu|/3) - 1 = mu**2*poly(mu**2)
(BODY_DB, degree 4, 2e-8 for m >= 1). Inverse tables of the far side, 8 mu
points by 16, uniform:
  x <= pi/2: bt_a(:,:,1) = x/xq (short way, mu in [0,1/pi]), xq the root of
  3xq**2/8 + 2*mu*xq = Q - 1, on xq/xq(pi/2); bt_a(:,:,2) = (x - x_b)/(xq - x_b)
  (long way, |mu| in [0,1/(2pi)]), xq = x_b + sqrt(8*(Q - Q_b)/3), on
  (xq - x_b)/(xq(pi/2) - x_b);
  x >= pi/2: bt_b = cos(x/2)/c0, c0 = Q**(-1/3), on c0/c0(pi/2), signed mu in
  [-1/(2pi),1/pi].
Every value tends to a finite limit at its ends, so the interpolated root
keeps its relative accuracy at x_b and at pi (within 3e-4 of min(x - x_b, pi - x)).

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


MR_M = [1, 2, 3, 4, 5, 6, 8, 12, 16, 32, 1000000]; MR_NT = 16; MR_NS = 12; MR_NC = 128
MR_TC = 1.45      # mr_q tau grid: tau = TM*sin(MR_TC*x)/sin(MR_TC), x uniform in [-1,1]


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
    tg = TM * np.sin(MR_TC * np.linspace(-1, 1, n)) / np.sin(MR_TC); lg = 2 * tg / (1 + tg * tg)
    sg = np.linspace(0, 1, MR_NS) ** 2                 # uniform in sqrt|sigma|
    tc = np.linspace(-TM, TM, MR_NC); lc = 2 * tc / (1 + tc * tc)
    V = np.zeros((n, MR_NS, 2, len(MR_M))); C = np.zeros((MR_NC, len(MR_M))); VS = np.zeros_like(C)
    for k, m in enumerate(MR_M):
        _, vg = tmin_v(lg, m)
        for h, br in enumerate((-1, 1)):
            S = br * np.clip(sg, 1e-7, 1 - 1e-9)
            L, SS = np.meshgrid(lg, S, indexing='ij')
            P = np.pi if br > 0 else 0.0
            V[:, :, h, k] = ((branch_v(L, SS, m) - vg[:, None]) / (np.abs(SS) * (P - vg[:, None])) - 1) / (1 - np.abs(SS))
        tm, vs = tmin_v(lc, m)
        # T scales as s_L**1.5 (canonical s_L = 2/(1+lam**2)).
        C[:, k] = np.cbrt(tm / m) / np.sqrt(2 / (1 + lc * lc)); VS[:, k] = vs
    return V, C, VS


TE0, TE1 = -0.9, 1 - 1e-6            # tau range of the elliptic table tabe
NE_T, NE_S = 40, 32                    # its tau and s points


def zc(tau):
    """Grid coordinate of tabe: sign(tau)*(1 - sqrt(1 - |tau|)), fine toward both ends."""
    return np.sign(tau) * (1 - np.sqrt(1 - np.abs(tau)))


def half_scales(lam):
    """tp, T_half = T(pi/2) and K = pi*y(pi)**1.5/sqrt2 for the canonical geometry."""
    a, b, y0, yend, tp = geom(lam)
    ypi = np.where(b >= 0, y0 + 2 * b, yend)
    yh = np.where(b >= 0, y0 + b, yend - b)
    return tp, np.sqrt(yh) * (yh * PI / (2 * S2) + a), PI * ypi ** 1.5 / S2


def ellipse(nt=NE_T, ns=NE_S):
    """Every zero-revolution ellipse with tau in [TE0, TE1], split at v = pi/2 (T_half):
    tabe(:,:,1) = (pi-v)/(s1*(K/T_half)**(1/3)) on s1 = (T_half/t)**(1/3) (pole side);
    tabe(:,:,2) = v**2/(s2*(s2+2*e0)) on s2 = (t-tp)/(T_half-tp), e0 = tp/(T_half-tp)
    (parabola side), which absorbs the sqrt(y) branch point's approach to the
    parabola as tau -> 1. Both stay finite at their ends. The tau grid is fine
    toward both ends; one table replaces the former tab (|tau| <= 0.9, split at
    s = 1/2) and tabp (tau >= 0.75) with 2,560 values instead of 4,096 and the
    same accuracy (0.35% of the singularity bound in the middle, 3.3% at
    tau = -0.9, 3.2e-5 toward tau = 1)."""
    z = np.linspace(zc(TE0), zc(TE1), nt); tau = np.sign(z) * (1 - (1 - np.abs(z)) ** 2)
    lam = 2 * tau / (1 + tau * tau)
    L, S = np.meshgrid(lam, np.linspace(0, 1, ns), indexing='ij')
    tp, th, K = half_scales(L)
    s1 = np.maximum(S, 1e-7)
    A = (PI - root_t(L, th / s1 ** 3)) / (s1 * np.cbrt(K / th))
    s2 = np.clip(S, 1e-7, 1 - 1e-12); e0 = tp / (th - tp)
    B = root_t(L, tp + s2 * (th - tp)) ** 2 / (s2 * (s2 + 2 * e0))
    return A, B


BB_DEG = 5                              # degree of each half of the body fit


def P_body(d):
    """T/|b|**1.5 near the long-way limit y_end -> 0, d = pi - v, without
    cancellation: in v = pi - d, (v - sin v cos v)/(4 sin(v/2)**3) - cos(v/2),
    with the series of v - sin(2v)/2 below v = 0.5."""
    v = PI - d; s = np.sin(v / 2)
    ser = np.zeros_like(v); term = np.ones_like(v)          # sum (-1)^(k+1) 4^k v^(2k+1)/(2(2k+1)!)
    for k in range(1, 12):
        term = term * (2 * v) ** 2 / ((2 * k) * (2 * k + 1)) if k > 1 else (2 * v) ** 3 / 6
        ser = ser + (-1) ** (k + 1) * term / 2
    num = np.where(v < .5, ser, v - np.sin(2 * v) / 2)
    return num / (4 * s ** 3) - np.cos(v / 2)


def long_body(deg=BB_DEG):
    """Long-way limit tau -> -1: d = pi - v from P_body(d) = t/|b|**1.5 in
    u = (t/|b|**1.5 - 1/3)/(pi/4 - 1/3) (P_body is 1/3 at the parabola and
    pi/4 at the pole). Least-squares polynomials of degree deg in the local
    variable of each half: v**2/u = sum bb(k,1) (4u-1)**k on [0,1/2] and
    d/(1-u) = sum bb(k,2) (4u-3)**k on [1/2,1]. Both quantities are smooth and
    finite at their ends; degree 5 is within 2.4e-5 of min(v,d)."""
    d = np.r_[np.linspace(1e-6, PI - 1e-6, 6000), PI - np.logspace(-6, -1, 400)]
    u = (P_body(d) - 1 / 3) / (PI / 4 - 1 / 3); v = PI - d
    lo = u <= .5
    P = np.polynomial.Polynomial
    c1 = P.fit(u[lo], v[lo] ** 2 / u[lo], deg, domain=[0, .5]).coef
    c2 = P.fit(u[~lo], d[~lo] / (1 - u[~lo]), deg, domain=[.5, 1]).coef
    return np.stack([c1, c2], axis=-1)


def root_t(lam, t):
    """Elliptic root v of T(v) = t by bisection."""
    a, b, y0, yend, tp = geom(lam)
    lo = np.zeros_like(lam); hi = np.full(lam.shape, PI)
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        f = T_v(mid, a, b, y0, yend) < t
        lo = np.where(f, mid, lo); hi = np.where(f, hi, mid)
    return 0.5 * (lo + hi)


BT_NM = 8; BT_NA = 16; MU_S = 1 / PI; MU_L = 1 / (2 * PI)


def body_Q(x, mu):
    return (1 + mu * (x + np.sin(x))) / np.cos(x / 2) ** 3


def body_inv(Qv, mu, lo, hi):
    """x on [lo, hi] (Q increasing there) with body_Q(x, mu) = Qv, by bisection."""
    lo = np.broadcast_to(np.asarray(lo, float), np.shape(Qv)).copy()
    hi = np.broadcast_to(np.asarray(hi, float), np.shape(Qv)).copy()
    for _ in range(100):
        m = .5 * (lo + hi); up = body_Q(m, mu) < Qv
        lo = np.where(up, m, lo); hi = np.where(up, hi, m)
    return .5 * (lo + hi)


def body_db_poly():
    """Coefficients (highest first) of (x_b/(8a/3) - 1)/a**2 in a**2, a = |mu|,
    x_b the long-way minimum: a = 1/(x + sin x + 4 cos(x/2)**3/(3 sin(x/2)))."""
    k = np.arange(40); a = MU_L * (1 - np.cos(PI * (k + .5) / 40)) / 2
    lo = np.full_like(a, 1e-12); hi = np.full_like(a, 2.0)
    for _ in range(200):
        d = .5 * (lo + hi); up = 1 / (d + np.sin(d) + 4 * np.cos(d / 2) ** 3 / (3 * np.sin(d / 2))) < a
        lo = np.where(up, d, lo); hi = np.where(up, hi, d)
    return np.polyfit(a * a, (.5 * (lo + hi) / (8 * a / 3) - 1) / a ** 2, 4)


def body_min(a, C):
    """x_b and Q(x_b) - 1 (cancellation-free), as the solver computes them."""
    db = 8 * a / 3 * (1 + a * a * np.polyval(C, a * a)); c = np.cos(db / 2)
    return db, (2 * np.sin(db / 4) ** 2 * (1 + c + c * c) - a * (db + np.sin(db))) / c ** 3


def body_tables():
    C = body_db_poly(); a = np.linspace(0, 1, BT_NA); qh = lambda mu: body_Q(PI / 2, mu)
    A = np.zeros((BT_NA, BT_NM, 2)); B = np.zeros((BT_NA, BT_NM))
    for k, mu in enumerate(np.linspace(0, MU_S, BT_NM)):          # short way
        eh = qh(mu) - 1; xh = 2 * eh / (2 * mu + np.sqrt(4 * mu * mu + 1.5 * eh)); xq = a * xh
        x = body_inv(1 + 3 * xq ** 2 / 8 + 2 * mu * xq, mu, 0, PI / 2 + 1e-9)
        A[:, k, 0] = np.where(a > 0, x / np.where(a > 0, xq, 1), 1)
    for k, am in enumerate(np.linspace(0, MU_L, BT_NM)):          # long way
        db, em = body_min(am, C); xh = db + np.sqrt(8 * (qh(-am) - 1 - em) / 3); xq = db + a * (xh - db)
        x = body_inv(1 + em + 3 * (xq - db) ** 2 / 8, -am, db, PI / 2 + 1e-9)
        c = np.cos(db / 2); s = np.sin(db / 2); n = 1 - am * (db + np.sin(db))
        L2 = am * np.sin(db) / n - (am * (1 + np.cos(db)) / n) ** 2 + .75 / c ** 2   # (log Q)''
        A[:, k, 1] = np.where(a > 0, (x - db) / np.where(a > 0, xq - db, 1), np.sqrt(.75 / ((1 + em) * L2)))
    for k, mu in enumerate(np.linspace(-MU_L, MU_S, BT_NM)):      # far side, both ways
        c0 = a * qh(mu) ** (-1 / 3)
        x = body_inv(np.where(a > 0, 1 / np.maximum(c0, 1e-300) ** 3, np.inf), mu, PI / 2 - 1e-9, PI)
        B[:, k] = np.where(a > 0, np.cos(x / 2) / np.maximum(c0, 1e-300), (1 + mu * PI) ** (1 / 3))
    return A, B, C


def lit(x):
    """Shortest literal for float32(x) (correctly rounded on input)."""
    x = np.float32(x)
    a = np.format_float_positional(x, unique=True, trim='.')
    b = np.format_float_scientific(x, unique=True, trim='.', exp_digits=1)
    r = min(a, b, key=len)
    return r.replace('0.', '.', 1) if r.startswith(('0.', '-0.')) else r


def arr(f, name, M):
    """Declare a protected array and initialize it with DATA statements of
    whole columns (at most 1000 values, under the standard's 255
    continuation lines). Values are stored in single precision (relative
    rounding 6e-8, far below every starter's error) as the shortest decimal
    that rounds back to the same float32."""
    shp = M.shape
    f.write(f"   real(real32),protected::{name}({','.join(str(x) for x in shp)})\n")
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
        f.write(", &\n".join("      " + ",".join(lit(x) for x in vals[k:k + 8]) for k in range(0, len(vals), 8)))
        f.write("/\n")


def write(out, n, EA, EB, LB, HS, HL, V, C, VS, BA, BB, BC):
    def block(f, name, mats):
        arr(f, name, np.stack(mats, axis=-1))
    with open(out, 'w') as f:
        f.write(f"   ! Generated by tools/gen_tables.py {n}; do not edit.\n")
        f.write(f"   ! Zero-revolution root tables on tau = lam/(1+sqrt(1-lam**2)).\n")
        f.write("   integer,private::ti_,tj_\n")
        f.write(f"   ! Ellipse, tau in [{TE0},1-{1-TE1:.0e}], grid uniform in z = sign(tau)*(1-sqrt(1-|tau|)),\n")
        f.write("   ! split at v=pi/2: tabe(:,:,1)=(pi-v)/(s1*(K/Th)**(1/3)), s1=(Th/t)**(1/3), Th=T(pi/2),\n")
        f.write("   ! K=pi*y(pi)**1.5/sqrt2; tabe(:,:,2)=v**2/(s2*(s2+2*e0)), s2=(t-tp)/(Th-tp), e0=tp/(Th-tp).\n")
        f.write(f"   integer,parameter::tabe_nt={EA.shape[0]},tabe_ns={EA.shape[1]}\n")
        f.write(f"   real(dp),parameter::tabe_z0={zc(TE0):.17e}_dp,tabe_z1={zc(TE1):.17e}_dp\n")
        block(f, "tabe", (EA, EB))
        f.write("   ! Long-way limit tau -> -1, body d=pi-v from T/|b|**1.5=th, u=(th-1/3)/(pi/4-1/3):\n")
        f.write("   ! v**2/u = sum bb(k,1)*(4u-1)**k on [0,1/2], d/(1-u) = sum bb(k,2)*(4u-3)**k on [1/2,1].\n")
        f.write(f"   integer,parameter::bb_n={LB.shape[0]}\n")
        f.write("   real(dp),protected::bb(bb_n,2)=reshape([" + ",".join(f"{x:.17e}_dp" for x in LB.T.ravel()) + "],[bb_n,2])\n")
        f.write(f"   ! Hyperbola, sqrt|tau| in [0,sqrt({TOP})], u = t/tp in [{UH},1]: tabh(:,:,1) = w**2/((1-u)*wn**2),\n")
        f.write("   ! wn**2 = L**2*(1-exp(-wb**2/L**2)), L = 1-log(u) (short way); tabh(:,:,2) = w**2/((1-u)*L**2) (long way).\n")
        f.write(f"   integer,parameter::tab_n={n}\n   real(dp),parameter::tab_tmax={TM}_dp\n")
        f.write(f"   real(dp),parameter::tabh_s1={np.sqrt(TOP):.17e}_dp,tabh_u0={UH}_dp\n")
        block(f, "tabh", (HS, HL))
        f.write(f"   ! Multiple revolutions (slices in mu = 1/m, last ~infinity): mr_q(tau,|sigma|,branch,k) = ((v-v*)/(|sigma|*(P-v*)) - 1)/(1-|sigma|),\n")
        f.write(f"   ! mr_c(tau,k) = cbrt(t_min/m)/sqrt(s_L), mr_s(tau,k) = v*; |tau| <= {TM}.\n")
        f.write(f"   integer,parameter::mr_nt={V.shape[0]},mr_ns={V.shape[1]},mr_nm={len(MR_M)},mr_nc={MR_NC}\n")
        f.write(f"   ! mr_q grid: tau = {TM}*sin(mr_tc*x)/sin(mr_tc), x uniform in [-1,1]; |sigma| = (j/(mr_ns-1))**2.\n")
        f.write(f"   real(dp),parameter::mr_tc={MR_TC}_dp,mr_tk={np.sin(MR_TC) / TM:.17e}_dp\n")
        f.write("   real(dp),parameter::mr_mu(mr_nm)=[" + ",".join(f"{1 / m:.17e}_dp" for m in MR_M) + "]\n")
        f.write("   integer,parameter::mr_m(mr_nm)=[" + ",".join(str(m) for m in MR_M) + "]\n")
        arr(f, "mr_q", V); arr(f, "mr_c", C); arr(f, "mr_s", VS)
        f.write("   ! Body-limit inverses (mu = 1/M short way, -1/M long way; see tools/gen_tables.py):\n")
        f.write("   ! bt_a(:,:,1) = x/xq, bt_a(:,:,2) = (x-x_b)/(xq-x_b) for x <= pi/2; bt_b = cos(x/2)/c0 beyond.\n")
        f.write(f"   integer,parameter::bt_na={BT_NA},bt_nm={BT_NM}\n")
        f.write("   real(dp),parameter::bt_db(5)=[" + ",".join(f"{x:.17e}_dp" for x in BC) + "]\n")
        arr(f, "bt_a", BA); arr(f, "bt_b", BB)

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    out = sys.argv[2] if len(sys.argv) > 2 else "src/excali_tables.inc"
    EA, EB = ellipse()
    LB = long_body()
    HS, HL = hyperbolic(n)
    V, C, VS = multi_rev(MR_NT)
    BA, BB, BC = body_tables()
    write(out, n, EA, EB, LB, HS, HL, V, C, VS, BA, BB, BC)
