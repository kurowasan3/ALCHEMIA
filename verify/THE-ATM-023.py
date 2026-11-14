# THE-ATM-023 検証：水素原子のエネルギー準位と発光・電離
from sympy import Rational as Q, N, oo
A = Q('2.18e-18'); hc = Q('1.99e-25'); NA = Q('6.02e23')
def dE(nu, nl): return A*(Q(1, nl**2) - (0 if nu is oo else Q(1, nu**2)))
def lam_nm(nu, nl): return hc/dE(nu, nl)*10**9
l42 = lam_nm(4, 2); print('4->2', N(l42, 5)); assert round(float(l42)) == 487
lines = {n: N(lam_nm(n, 2), 5) for n in range(3, 9)}; print(lines)
vis = [n for n, l in lines.items() if 400 <= l <= 700]; print('visible', vis); assert vis == [3, 4, 5, 6]
linf = lam_nm(oo, 2); print('series limit', N(linf, 5)); assert round(float(linf)) == 365
I1 = A*NA/1000; I2 = A/4*NA/1000
print('I(n=1) kJ/mol', N(I1, 5), ' I(n=2)', N(I2, 5))
assert round(float(I1), -1) == 1310 and round(float(I2)) == 328
Eph = hc/Q('300e-9'); print('photon 300 nm J', N(Eph, 5), ' vs ', N(A/4), N(A))
assert A/4 < Eph < A
# ライマン系列の最長波長（2->1）
print('Lyman longest', N(hc/dE(2, 1)*1e9, 5))
print('OK')
