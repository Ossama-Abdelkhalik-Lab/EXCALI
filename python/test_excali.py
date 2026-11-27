"""pytest for python/excali.py. Needs the C library: set EXCALI_LIBRARY or
build into ../build. Compares the Python wrapper with the C library called
directly through ctypes (bit-identical), the batch with per-case solves (to
rounding), checks no-copy passing and concurrent calls from threads."""
import ctypes
import math
import os
import sys
import threading

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import excali  # noqa: E402

MU = 398600.4418


def cases(n, seed=1):
    rng = np.random.default_rng(seed)
    q = rng.normal(size=(n, 4))
    q /= np.linalg.norm(q, axis=1)[:, None]
    th = rng.uniform(0.05, 2 * math.pi - 0.05, n)
    rr1, rr2 = rng.uniform(6500, 42000, n), rng.uniform(6500, 42000, n)
    r1 = np.zeros((n, 3)); r2 = np.zeros((n, 3))
    for i in range(n):
        a, b, c, d = q[i]
        R = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                      [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                      [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])
        r1[i] = R @ [rr1[i], 0, 0]
        r2[i] = R @ [rr2[i]*math.cos(th[i]), rr2[i]*math.sin(th[i]), 0]
    s = 0.5 * (rr1 + rr2 + np.linalg.norm(r2 - r1, axis=1))
    unit = np.sqrt(s**3 / MU)
    nrev = np.where(rng.uniform(size=n) < 0.3, rng.integers(1, 6, n), 0).astype(np.intc)
    tof = np.where(nrev == 0, unit * 10**rng.uniform(-1.5, 1, n),
                   unit * nrev * math.pi / math.sqrt(2) * rng.uniform(0.9, 4, n))
    direction = np.where(rng.uniform(size=n) < 0.5, 1, -1).astype(np.intc)
    branch = np.where(rng.uniform(size=n) < 0.5, 1, -1).astype(np.intc)
    return r1, r2, tof, direction, nrev, branch


def test_solve_matches_c():
    lib = excali.load()
    r1, r2, tof, d, m, br = cases(300)
    nsolved = 0
    for i in range(300):
        v1, v2, st = excali.solve(r1[i], r2[i], tof[i], MU, d[i], m[i], br[i])
        w1, w2, s = np.empty(3), np.empty(3), ctypes.c_int()
        a, b = r1[i].copy(), r2[i].copy()
        lib.excali_solve(a.ctypes.data, b.ctypes.data, tof[i], MU, int(d[i]), int(m[i]), int(br[i]), w1.ctypes.data, w2.ctypes.data,
                         ctypes.byref(s))
        assert st == s.value
        np.testing.assert_array_equal(v1, w1)
        np.testing.assert_array_equal(v2, w2)
        nsolved += st == 0
    assert nsolved > 200


def test_batch_matches_solve():
    r1, r2, tof, d, m, br = cases(1000, seed=2)
    v1, v2, st, half = excali.batch(r1, r2, tof, MU, d, m, br)
    for i in range(1000):
        w1, w2, s = excali.solve(r1[i], r2[i], tof[i], MU, d[i], m[i], br[i])
        assert s == st[i]
        if s in (0, 1, 2):
            assert np.max(np.abs(np.r_[v1[i] - w1, v2[i] - w2])) <= 1e-12 * np.linalg.norm(w1)
        else:
            assert np.isnan(v1[i]).all()


def test_batch_mu1_and_scalars():
    r1, r2, tof, d, m, br = cases(64, seed=3)
    t1 = tof * math.sqrt(MU)  # mu = 1 time units
    a = excali.batch(r1, r2, t1, 1.0, d, m, br)
    b = excali.batch(r1, r2, tof, MU, d, m, br)
    np.testing.assert_array_equal(a[2], b[2])
    ok = a[2] == 0
    np.testing.assert_allclose(a[0][ok] * math.sqrt(MU), b[0][ok], rtol=1e-12)
    c = excali.batch(r1, r2, t1, direction=1, nrev=0)   # broadcast scalars
    assert c[0].shape == (64, 3) and c[2].dtype == np.intc


def test_zero_copy():
    r1, r2, tof, d, m, br = cases(100, seed=4)
    assert r1.flags.c_contiguous and r1.dtype == np.float64
    assert excali.data_address(r1) == r1.ctypes.data          # Fortran sees the same memory
    out = (np.empty((100, 3)), np.empty((100, 3)), np.empty(100, np.intc), np.empty(100, np.intc))
    res = excali.batch(r1, r2, tof, MU, d, m, br, out=out)
    assert all(x is y for x, y in zip(res, out))                # outputs written in place
    assert excali._rows(r1, 100, "r1") is r1                    # no conversion copy on input
    assert excali._ints(d, 100, "d") is d
    # A Fortran-ordered or float32 input is converted (a copy), with identical results.
    f = np.asfortranarray(r1)
    res2 = excali.batch(f, r2, tof, MU, d, m, br)
    np.testing.assert_array_equal(res2[0], res[0])


def test_invalid_inputs():
    v1, v2, st = excali.solve([1, 0, 0], [0, 1, 0], -1.0)
    assert st == -1 and np.isnan(v1).all()
    v1, v2, st = excali.solve([1, 0, 0], [-1, 0, 0], 3.0)
    assert st == -2
    v1, v2, st = excali.solve([1, 0, 0], [-1, 0, 0], 3.0, plane_normal=[0, 0, 1])
    assert st == 0
    r1, r2, tof, d, m, br = cases(5, seed=5)
    v1, v2, st, h = excali.batch(r1, r2, tof, -1.0, d, m, br)
    assert (st == -1).all()
    with pytest.raises(ValueError):
        excali.batch(r1[:, :2], r2, tof)
    empty = excali.batch(np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0))
    assert empty[0].shape == (0, 3)


def test_threads():
    r1, r2, tof, d, m, br = cases(20000, seed=6)
    ref = excali.batch(r1, r2, tof, MU, d, m, br)
    results = [None] * 4

    def work(k):
        results[k] = excali.batch(r1, r2, tof, MU, d, m, br)

    th = [threading.Thread(target=work, args=(k,)) for k in range(4)]
    for t in th:
        t.start()
    for t in th:
        t.join()
    for r in results:
        np.testing.assert_array_equal(r[2], ref[2])
        np.testing.assert_array_equal(r[0], ref[0])
