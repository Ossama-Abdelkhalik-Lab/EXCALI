from fractions import Fraction as Fr
from math import factorial
N=24
# v - sin v cos v = v - sin(2v)/2 ; divide by v^3 -> series in X
num=[]
for n in range(1,N+2):  # term (2v)^(2n+1)/(2(2n+1)!) * (-1)^(n+1), power v^(2n+1)
    num.append(Fr((-1)**(n+1)*2**(2*n+1),2*factorial(2*n+1)))  # coefficient of v^(2n+1); /v^3 -> X^(n-1)
# sin v / v = sum (-1)^n X^n/(2n+1)!
s=[Fr((-1)**n,factorial(2*n+1)) for n in range(N+1)]
def mul(a,b):
    return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(N+1)]
s3=mul(mul(s,s),s)
# num/s3
q=[]
for k in range(N+1):
    q.append((num[k]-sum(q[i]*s3[k-i] for i in range(k)))/s3[0])
print([str(x) for x in q[:10]])
for x in q[:22]: print('%.17e'%float(x), float(q[q.index(x)+1]/x) if q.index(x)<21 else '')
print()
M=22
co=[q[k] for k in range(M)][::-1]
s=",".join("%.17e_dp"%float(x) for x in co).replace("e-0","e-").replace("e+00","")
items=["%.16e_dp"%float(x) for x in co]
lines=[]; cur=[]
for it in items:
    cur.append(it)
    if len(cur)==4: lines.append(",".join(cur)); cur=[]
if cur: lines.append(",".join(cur))
print("   real(dp),parameter :: pall(%d)=[ &\n      "%M + ", &\n      ".join(lines)+"]")
