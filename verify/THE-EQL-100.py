# THE-EQL-100 検証：光路中の NO2 量の比から N2O4 ⇄ 2NO2 の Kp を逆算し，比のとりうる範囲を求める
import sympy as sp
r, y = sp.symbols('r y', positive=True)
PA, PB = sp.Rational(1, 10), sp.Rational(1, 20)
# 光路中 NO2 ∝ y·P·L，P·L が等しいので比 = yB/yA
eq = sp.Eq(PA*y**2/(1-y), PB*(r*y)**2/(1-r*y))
ysol = sp.solve(eq, y); ysol = [s for s in ysol if s != 0]
print('y_A(r) =', sp.simplify(ysol[0]))
yf = sp.simplify(ysol[0])
assert sp.simplify(yf - (2 - r**2)/(r*(2 - r))) == 0
yv = yf.subs(r, sp.Rational(1236, 1000)); print('y_A =', float(yv)); assert round(float(yv), 2) == 0.50
Kp = PA*yv**2/(1-yv); print('Kp =', float(Kp)); assert f"{float(Kp):.2g}" == '0.05'
# 範囲：y→0 で r→√2，y→1 で r→1
print('limits', sp.solve(sp.Eq(2 - r**2, 0), r), sp.solve(sp.Eq(2 - r**2, r*(2 - r)), r))
# r = 1.50 は y<0 となり不可能
print('r=1.5 ->', float(yf.subs(r, sp.Rational(3, 2))))
assert yf.subs(r, sp.Rational(3, 2)) < 0
# 予測：P = 0.025 atm, L = 40 cm
yc = [s for s in sp.solve(sp.Eq(y**2/(1-y), sp.Rational(5, 100)/sp.Rational(25, 1000)), y) if 0 < s < 1][0]
print('y_C =', float(yc), ' ratio to A', float(yc/sp.Rational(1, 2)))
assert round(float(yc/sp.Rational(1, 2)), 2) == 1.46
