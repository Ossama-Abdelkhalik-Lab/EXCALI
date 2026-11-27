"""Python interface to EXCALI through its C library (include/excali.h).

ctypes and NumPy only. Arrays of 3-vectors are (n, 3) C-contiguous float64,
the same bytes as the Fortran (3, n) arrays, so batches are passed without a
copy; integer arrays are C int (int32 on every supported platform). ctypes
releases the GIL during each call.

    import excali
    v1, v2, status = excali.solve(r1, r2, tof, mu=398600.4418)
    v1, v2, status, half = excali.batch(r1s, r2s, tofs, mu=398600.4418)

The library is found from the EXCALI_LIBRARY environment variable, or as
libexcali_c.{dylib,so,dll} in the given build directory or ../build.
"""
import ctypes
import os
import sys

import numpy as np

STATUS = {0: "solved", 1: "coalesced roots at minimum time", 2: "solved, low confidence",
          5: "time below multi-revolution minimum", -1: "invalid input",
          -2: "degenerate geometry or missing antipodal plane", -3: "nonunique collinear geometry",
          -4: "numerical failure"}

_lib = None


def load(path=None):
    """Load the C library (once). path: the library file or its directory."""
    global _lib
    if _lib is not None and path is None:
        return _lib
    name = {"darwin": "libexcali_c.dylib", "win32": "excali_c.dll"}.get(sys.platform, "libexcali_c.so")
    cands = []
    if path:
        cands.append(path if os.path.isfile(path) else os.path.join(path, name))
    if os.environ.get("EXCALI_LIBRARY"):
        cands.append(os.environ["EXCALI_LIBRARY"])
    cands.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build", name))
    for c in cands:
        if os.path.isfile(c):
            lib = ctypes.CDLL(os.path.abspath(c))
            break
    else:
        raise OSError(f"EXCALI C library not found; tried {cands}. Build it (cmake) or set EXCALI_LIBRARY.")
    # Raw addresses (c_void_p): ndpointer's per-argument checks cost
    # microseconds per call; the wrappers below validate shapes and dtypes.
    vp, ip = ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)
    c_d, c_i = ctypes.c_double, ctypes.c_int
    lib.excali_solve.argtypes = [vp, vp, c_d, c_d, c_i, c_i, c_i, vp, vp, ip]
    lib.excali_solve_plane.argtypes = [vp, vp, c_d, c_d, c_i, c_i, c_i, vp, vp, vp, ip]
    lib.excali_batch.argtypes = [c_i] + [vp] * 10
    lib.excali_batch_mu.argtypes = [c_i, vp, vp, vp, c_d] + [vp] * 7
    lib.excali_c_data_address.argtypes = [vp]
    lib.excali_c_data_address.restype = ctypes.c_ssize_t
    for f in ("excali_solve", "excali_solve_plane", "excali_batch", "excali_batch_mu"):
        getattr(lib, f).restype = None
    _lib = lib
    return lib


def _vec(x, name):
    a = x if (type(x) is np.ndarray and x.dtype == np.float64 and x.flags.c_contiguous) \
        else np.ascontiguousarray(x, dtype=np.float64)
    if a.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {a.shape}")
    return a


def solve(r1, r2, tof, mu=1.0, direction=1, nrev=0, branch=-1, plane_normal=None):
    """One transfer. Returns (v1, v2, status); velocities are NaN unless
    status is 0, 1 or 2."""
    lib = load()
    v = np.empty(6)
    st = ctypes.c_int()
    a, b = _vec(r1, "r1"), _vec(r2, "r2")
    p = v.ctypes.data
    args = (a.ctypes.data, b.ctypes.data, float(tof), float(mu), int(direction), int(nrev), int(branch))
    if plane_normal is None:
        lib.excali_solve(*args, p, p + 24, ctypes.byref(st))
    else:
        lib.excali_solve_plane(*args, _vec(plane_normal, "plane_normal").ctypes.data, p, p + 24, ctypes.byref(st))
    return v[:3], v[3:], st.value


def _rows(x, n, name):
    """(n,3) float64 C-contiguous: the input itself when it already is one."""
    a = np.asarray(x)
    if a.dtype != np.float64 or not a.flags.c_contiguous:
        a = np.ascontiguousarray(a, dtype=np.float64)
    if a.shape != (n, 3):
        raise ValueError(f"{name} must have shape ({n}, 3), got {a.shape}")
    return a


def _ints(x, n, name):
    a = np.asarray(x)
    if a.ndim == 0:
        return np.full(n, int(a), dtype=np.intc)
    if a.dtype != np.intc or not a.flags.c_contiguous:
        a = np.ascontiguousarray(a, dtype=np.intc)
    if a.shape != (n,):
        raise ValueError(f"{name} must be a scalar or have shape ({n},), got {a.shape}")
    return a


def batch(r1, r2, tof, mu=1.0, direction=1, nrev=0, branch=-1, out=None):
    """n transfers sharing one mu. r1, r2: (n,3); tof: (n,); direction, nrev,
    branch: scalars or (n,). Float64 C-contiguous inputs and int32 index
    arrays are passed without a copy. out: optional (v1, v2, status, half)
    arrays to fill. Returns (v1 (n,3), v2 (n,3), status (n,), half (n,))."""
    lib = load()
    t = np.asarray(tof)
    if t.dtype != np.float64 or not t.flags.c_contiguous:
        t = np.ascontiguousarray(t, dtype=np.float64)
    if t.ndim != 1:
        raise ValueError("tof must be one-dimensional")
    n = t.shape[0]
    a, b = _rows(r1, n, "r1"), _rows(r2, n, "r2")
    d, m, br = _ints(direction, n, "direction"), _ints(nrev, n, "nrev"), _ints(branch, n, "branch")
    if out is None:
        out = (np.empty((n, 3)), np.empty((n, 3)), np.empty(n, dtype=np.intc), np.empty(n, dtype=np.intc))
    v1, v2, st, hf = out
    if (v1.shape != (n, 3) or v2.shape != (n, 3) or st.shape != (n,) or hf.shape != (n,)
            or v1.dtype != np.float64 or v2.dtype != np.float64 or st.dtype != np.intc or hf.dtype != np.intc
            or not all(x.flags.c_contiguous and x.flags.writeable for x in out)):
        raise ValueError("out must be (v1 (n,3) float64, v2 (n,3) float64, status (n,) int32, half (n,) int32), "
                         "C-contiguous and writeable")
    if n == 0:
        return out
    ptr = [x.ctypes.data for x in (a, b, t, d, m, br, v1, v2, st, hf)]
    if float(mu) == 1.0:
        lib.excali_batch(n, *ptr)
    else:
        lib.excali_batch_mu(n, *ptr[:3], float(mu), *ptr[3:])
    return out


def data_address(a):
    """The address the Fortran library receives for array a (a no-copy check)."""
    return load().excali_c_data_address(a.ctypes.data)
