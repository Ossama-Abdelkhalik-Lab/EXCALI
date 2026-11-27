"""Compare host_driver output (%a) with fortran_driver output (hex bits): statuses and ulp distances."""
import struct, sys
import numpy as np
def load_c(f):
    d = {}
    for l in open(f):
        p = l.split(); d[int(p[0])] = (int(p[1]), [float.fromhex(x) for x in p[2:]])
    return d
def load_f(f):
    d = {}
    for l in open(f):
        p = l.split(); d[int(p[0])] = (int(p[1]), [struct.unpack('>d', bytes.fromhex(x))[0] for x in p[2:]])
    return d
def ulps(a, b):
    ia = np.array(a).view(np.int64); ib = np.array(b).view(np.int64)
    ia = np.where(ia < 0, np.int64(-2**63) - ia, ia); ib = np.where(ib < 0, np.int64(-2**63) - ib, ib)
    return np.abs(ia - ib)
c, f = load_c(sys.argv[1]), load_f(sys.argv[2])
assert c.keys() == f.keys(), "row sets differ"
st_diff = [k for k in c if c[k][0] != f[k][0]]
u, rel = [], []
for k in c:
    if c[k][0] in (0, 2) and f[k][0] in (0, 2):
        a = np.array(c[k][1]); b = np.array(f[k][1])
        u.append(ulps(a, b).max())
        rel.append(max(np.linalg.norm(a[:3] - b[:3]), np.linalg.norm(a[3:] - b[3:])) / np.linalg.norm(b[:3]))
u = np.array(u); rel = np.array(rel)
print(f"{sys.argv[1]}: rows {len(c)}, status mismatches {len(st_diff)} {st_diff[:5]}, solved {len(u)}")
if len(u):
    print(f"  max-ulp per case: identical {np.mean(u == 0)*100:.2f}%  <=4 ulp {np.mean(u <= 4)*100:.2f}%  p99 {np.quantile(u, .99):.0f}  max {u.max()}")
    print(f"  relative velocity difference: p50 {np.median(rel):.2e} p99 {np.quantile(rel, .99):.2e} max {rel.max():.2e}")
