# THE-SOL-068 検証：水洗による水素の精製（全圧一定でのヘンリーの法則）
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, Rational as Q, nsimplify
assert check_equation('CH4 + 2H2O -> CO2 + 4H2') is None
nCO2, nH2 = Q(1, 10), Q(4, 10)       # CH4 0.10 mol から
s = Q(75, 1000)                      # CO2 分圧 1.0e5 Pa で水 1 L に 0.075 mol
def purity(Vw, P):                   # P: 全圧（1.0e5 Pa 単位）
    y = symbols('y', positive=True)  # 気相中の CO2 のモル分率
    sol = [r for r in solve(y*(nCO2 - s*Vw*P*y + nH2) - (nCO2 - s*Vw*P*y), y) if 0 < r < 1]
    assert len(sol) == 1
    return 1 - sol[0]
p1 = purity(5, 1); print('(1) 5.0 L, 1.0e5 Pa: purity =', float(p1))
assert round(float(p1)*100) == 88
# (2) 95 % にする水の量
y = Q(5, 100); gasCO2 = y/(1-y)*nH2; Vw = (nCO2 - gasCO2)/(s*y)
print('(2) water =', float(Vw)); assert round(float(Vw)) == 21
# (3) 全圧 2 倍・水 5 L と 全圧 1 倍・水 10 L
pa, pb = purity(5, 2), purity(10, 1)
print('(3)', float(pa), float(pb)); assert pa == pb and round(float(pa)*100) == 92
# 溶けた CO2 の量と気相の CO2
print('dissolved CO2 at (1):', float(s*5*(1-p1)))
print('OK')
