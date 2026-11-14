# THE-EQL-115 検証：2価/1価イオン交換の平衡と希釈の効果
from sympy import symbols, solve, sqrt, nsimplify
K, Q = 1.0, 2.0
r = symbols('r', positive=True)
def ca_frac(Na, Ca):
    a = K*Ca/Na**2               # s = a r^2
    sol = [float(x) for x in solve(r + 2*a*r**2 - Q, r) if 0 < float(x) <= Q]
    rr = sol[0]; s = a*rr**2
    return 2*s/Q, rr, s
f1 = ca_frac(0.010, 0.0010); f2 = ca_frac(0.0010, 0.00010); f3 = ca_frac(2.0, 0.0010)
print(f1, f2, f3)
assert round(f1[0]*100) == 85 and round(f2[0]*100) == 95 and round(f3[0]*100, 2) == 0.10
# 1:1 交換では希釈で比が変わらない
Kp = 2.5
for dil in (1, 10, 100):
    print('1:1 ratio', Kp*(0.0010/dil)/(0.010/dil))
# (5) Ca の当量分率 50 % となる [Na+]
Na50 = (K*0.0010/(0.5/1.0**2))**0.5; print('Na for 50%', Na50); assert round(Na50, 3) == 0.045
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
assert check_equation('CaCl2 -> Ca^2+ + 2Cl-') is None
print('OK')
