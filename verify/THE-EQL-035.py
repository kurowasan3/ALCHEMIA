# THE-EQL-035 検証：黒鉛と水蒸気の連立平衡（C + H2O <=> CO + H2，CO + H2O <=> CO2 + H2）
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, symbols, solve, nsolve, N, sqrt
for eq in ['C + H2O -> CO + H2', 'CO + H2O -> CO2 + H2', 'C + 2H2O -> CO2 + 2H2']:
    assert check_equation(eq) is None, eq
# 観測：H2O 1.00 mol の 60 % が消費，CO:CO2 = 3:1，全圧 2.0e5 Pa（一定）
x, y = symbols('x y')
s = solve([x + y - Q(6,10), (x - y) - 3*y], [x, y]); x0, y0 = s[x], s[y]
n = dict(H2O=1-x0-y0, CO=x0-y0, CO2=y0, H2=x0+y0); nt = 1 + x0
print(n, nt)
P = 2*10**5
Kp2 = n['CO2']*n['H2']/(n['CO']*n['H2O'])
Kp1 = n['CO']*n['H2']/n['H2O']*P/nt
print('Kp1', N(Kp1), 'Kp2', Kp2)
assert Kp2 == Q(1,2) and abs(N(Kp1) - 7.3e4) < 0.05e4
# (2) 転化率 80 % となる全圧 P'
s2 = solve([x + y - Q(8,10), y*(x+y) - Kp2*(x-y)*(1-x-y)], [x, y], dict=True)
s2 = [t for t in s2 if 0 < t[y] < t[x]]; x1, y1 = s2[0][x], s2[0][y]
P2 = Kp1 / ((x1-y1)*(x1+y1)/(1-x1-y1)/(1+x1))
print('x1,y1', x1, y1, "P'", N(P2), 'CO:CO2', (x1-y1)/y1)
assert (x1, y1) == (Q(72,100), Q(8,100)) and abs(N(P2) - 4.9e4) < 0.05e4 and (x1-y1)/y1 == 8
# 解の一意性（P' を与えたときの連立方程式の正の解が1つ）を数値で確認
X, Y = symbols('X Y', positive=True)
eqs = [ (X-Y)*(X+Y)/(1-X-Y)*P2/(1+X) - Kp1, Y*(X+Y) - Kp2*(X-Y)*(1-X-Y) ]
roots = set()
for gx in [0.3, 0.5, 0.7, 0.75, 0.85]:
    for gy in [0.02, 0.05, 0.1, 0.2]:
        try:
            r = nsolve(eqs, [X, Y], [gx, gy])
            if 0 < r[1] < r[0] and r[0] + r[1] < 1: roots.add((round(float(r[0]),6), round(float(r[1]),6)))
        except Exception: pass
print('roots', roots); assert roots == {(0.72, 0.08)}
# (3) 黒鉛 0.50 mol しかないとき（P' で）
assert x1 > Q(1,2)           # 黒鉛が不足する
yy = symbols('yy', positive=True)
ys = [r for r in solve(yy*(Q(1,2)+yy) - Kp2*(Q(1,2)-yy)**2, yy) if 0 < r < Q(1,2)][0]
print('y', ys, N(ys))
assert ys == -1 + sqrt(Q(5,4)) and abs(N(ys) - 0.118) < 0.001
nt3 = Q(3,2)
Q1 = (Q(1,2)-ys)*(Q(1,2)+ys)/(Q(1,2)-ys)*P2/nt3
print('Q1', N(Q1), '< Kp1', N(Kp1)); assert Q1 < Kp1
print('組成', N(Q(1,2)-ys,3), N(Q(1,2)-ys,3), N(ys,3), N(Q(1,2)+ys,3))
print('OK')
