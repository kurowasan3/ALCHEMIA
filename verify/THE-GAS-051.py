# THE-GAS-051 検証：2種の単分子膜で仕切り板をはさんだときのつり合い（2次元の状態方程式）
import sympy as sp
R = sp.Rational
c = R(40, 10)*sp.Integer(10)**-21     # J
W, L = R(20, 100), R(60, 100)        # m
N1, N2 = R(12, 10)*sp.Integer(10)**17, R(10, 10)*sp.Integer(10)**17
b1, b2 = R(20, 10)*sp.Integer(10)**-19, R(25, 10)*sp.Integer(10)**-19
P = lambda A, b: c/(A - b)
x0 = L/2
A1, A2 = W*x0/N1, W*(L-x0)/N2
print('initial A', float(A1), float(A2), ' P', float(P(A1, b1)), float(P(A2, b2)))
assert P(A1, b1) > P(A2, b2)
x = sp.symbols('x', positive=True)
xs = sp.solve(sp.Eq(P(W*x/N1, b1), P(W*(L-x)/N2, b2)), x)
xs = [s for s in xs if 0 < s < L]; assert len(xs) == 1; xe = xs[0]
print('x_eq =', float(xe), ' shift (cm) =', float((xe - x0)*100))
assert round(float(xe), 3) == 0.314 and round(float((xe-x0)*100), 1) == 1.4
Pe = P(W*xe/N1, b1); print('P_eq =', float(Pe)); assert f"{float(Pe):.2g}" == '0.012'
# 排除面積なし（b=0）の場合
xi = sp.solve(sp.Eq(c/(W*x/N1), c/(W*(L-x)/N2)), x)[0]
print('ideal x =', float(xi), N1/(N1+N2)*L); assert round(float(xi), 3) == 0.327
Pi = c/(W*xi/N1); print('ideal P', float(Pi))
# 線形方程式による手計算の形
lhs = W*x/N1 - b1 - (W*(L-x)/N2 - b2)
print(sp.solve(lhs, x))
