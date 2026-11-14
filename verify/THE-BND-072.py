# THE-BND-072 検証：閃亜鉛鉱型とウルツ鉱型の 1 組あたりの体積
from sympy import sqrt, Rational as Q, symbols, simplify, N
d = symbols('d', positive=True)
v_zb = (sqrt(2)*d)**3/4
v_wz = (sqrt(3)/2*d**2)*sqrt(Q(8, 3))*d/2
assert simplify(v_zb - d**3/sqrt(2)) == 0 and simplify(v_wz - d**3/sqrt(2)) == 0
v = float((d**3/sqrt(2)).subs(d, Q(41, 10**9)))
rho = 144/(6.0e23*v); print(v, rho); assert round(rho, 1) == 4.9
# 問題で与えた近似値でも
assert round(144/(6.0e23*0.0689e-21/1.41), 1) == 4.9
print('OK')
