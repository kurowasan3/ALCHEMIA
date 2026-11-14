# POL-SUG-030 検証：スクロースの酸加水分解（転化）を旋光度で追跡する一次反応
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, log, nsimplify, symbols, solve, N
assert check_equation('C12H22O11 + H2O -> C6H12O6 + C6H12O6') is None
Msuc = 12*12 + 22*1 + 11*16; Mmono = 6*12 + 12 + 6*16
assert (Msuc, Mmono) == (342, 180)
aS, aG, aF, aM = R(665, 10), R(527, 10), R(-924, 10), R(138)
l = 2
f = R(Mmono, Msuc)  # スクロース 1 g から各単糖 180/342 g
def alpha(cs, cm, x):  # x: スクロースの反応した割合
    return l*(aS*cs*(1 - x) + (aG + aF)*f*cs*x + aM*cm)
# 実験1
c0 = R(2, 10)
a0, ainf = alpha(c0, 0, 0), alpha(c0, 0, 1)
print('a0', N(a0), 'ainf', N(ainf, 6)); assert a0 == R(266, 10) and round(float(ainf), 2) == -8.36
x0 = solve(alpha(c0, 0, symbols('x')), symbols('x'))[0]
print('x at 0', N(x0, 6)); assert round(float(x0)*100, 1) == 76.1
# (3) 30 分で +9.12° → 残存率 1/2
frac = (R(912, 100) - ainf)/(a0 - ainf); print('frac30', N(frac, 6)); assert abs(float(frac) - 0.5) < 2e-4
k = log(2)/30; print('k', N(k, 4)); assert round(float(k), 4) == 0.0231
# (4) 60 分：残存率 1/4
a60 = alpha(c0, 0, R(3, 4)); print('a60', N(a60, 5)); assert round(float(a60), 2) == 0.38 and a60 > 0
t0 = N(log(1/(1 - x0))/k, 5); print('t(alpha=0)', t0); assert 60 < t0 < 65
# (7) 初濃度を 2 倍：x0 は c0 によらない
assert solve(alpha(2*c0, 0, symbols('x')), symbols('x'))[0] == x0
# 実験2
cs, cm = symbols('cs cm')
sol = solve([alpha(cs, cm, 0) - R(271, 10), alpha(cs, cm, 1) - R(962, 100)], [cs, cm])
print('cs,cm', N(sol[cs], 5), N(sol[cm], 5)); assert round(float(sol[cs]), 3) == 0.100 and round(float(sol[cm]), 4) == 0.0500
b0, binf = alpha(R(1, 10), R(5, 100), 0), alpha(R(1, 10), R(5, 100), 1)
b30 = alpha(R(1, 10), R(5, 100), R(1, 2)); print('b0,binf,b30', N(b0, 4), N(binf, 4), N(b30, 4))
assert round(float(b0), 1) == 27.1 and round(float(binf), 2) == 9.62 and round(float(b30), 1) == 18.4
# 実験2 の残存率の式はマルトースによらない：(b30-binf)/(b0-binf) = 1/2
assert (b30 - binf)/(b0 - binf) == R(1, 2)
assert binf > 0  # 旋光度は 0 にならない
print('POL-SUG-030 OK')
