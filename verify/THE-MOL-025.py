# THE-MOL-025 検証：シュウ酸カルシウム一水和物とシュウ酸マグネシウム二水和物の混合物の熱重量分析
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, symbols, solve
eqs = ['CaC2O4.H2O -> CaC2O4 + H2O', 'MgC2O4.2H2O -> MgC2O4 + 2H2O', 'CaC2O4 -> CaCO3 + CO',
       'MgC2O4 -> MgO + CO + CO2', 'CaCO3 -> CaO + CO2', 'MgC2O4 -> MgCO3 + CO', 'MgCO3 -> MgO + CO2', '2CO + O2 -> 2CO2']
for e in eqs:
    r = check_equation(e); assert r is None, (e, r)
H, C, O, Mg, Ca = 1, 12, 16, 24, 40
M_Ca = Ca + 2*C + 4*O + 2*H + O          # CaC2O4.H2O
M_Mg = Mg + 2*C + 4*O + 2*(2*H + O)      # MgC2O4.2H2O
assert (M_Ca, M_Mg) == (146, 148)
a0, b0 = Q(4,10), Q(3,10)    # mmol
m0 = M_Ca*a0 + M_Mg*b0
L1 = 18*a0 + 36*b0; L2 = 28*a0 + 72*b0; L3 = 44*a0
res = 56*a0 + 40*b0
print('m0', m0, 'losses', L1, L2, L3, 'residue', res)
assert m0 - L1 - L2 - L3 == res
assert (float(m0), float(L1), float(L2), float(L3)) == (102.8, 18.0, 32.8, 17.6)
# 仮説X：Mg も第2段階で MgCO3，第3段階で MgO
a, b = symbols('a b')
sX = solve([18*a + 36*b - L1, 44*(a+b) - L3], [a, b]); print('仮説X', sX)
assert sX[a] < 0
# 仮説Y（正しい）：第3段階は CaCO3 のみ
sY = solve([44*a - L3, 18*a + 36*b - L1], [a, b]); print('仮説Y', sY)
assert sY == {a: a0, b: b0}
assert 28*sY[a] + 72*sY[b] == L2 and M_Ca*sY[a] + M_Mg*sY[b] == m0
# 質量%と第2段階の気体の体積比
pct = M_Ca*a0/m0*100; print('CaC2O4.H2O %', float(pct))
assert abs(float(pct) - 56.8) < 0.05
print('CO:CO2 =', (a0+b0), ':', b0)   # 7:3
# 各プラトーの質量
print('plateaus', float(m0-L1), float(m0-L1-L2), float(res))
# 仮説X で第2段階の減少量を予測すると（a+b)*28
print('OK')
