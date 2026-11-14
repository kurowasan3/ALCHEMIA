# THE-RDX-036 検証：ヨウ素滴定による YBa2Cu3Oy の酸素量（銅の平均酸化数）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
import sympy as sp
y = sp.symbols('y')
M = sp.Rational('88.9') + 2*sp.Rational('137.3') + 3*sp.Rational('63.5') + 16*y
print('M =', M)
# Cu 1 個あたり Na2S2O3 は (1+p) mol，3(2+p) = 2y - 7 → 式量あたり 3+3p = 2y - 10
p = sp.symbols('p')
S_per = sp.expand(3*(1+p)).subs(p, sp.solve(sp.Eq(3*(2+p) + 3 + 4, 2*y), p)[0])
print('S2O3 per formula =', S_per); assert sp.simplify(S_per - (2*y - 10)) == 0
W = sp.Rational(1, 10); C = sp.Rational(2, 100); V = sp.Rational(286, 10)/1000
sol = sp.solve(sp.Eq(C*V, W*(2*y-10)/M), y)[0]
print('y =', float(sol)); assert round(float(sol), 2) == 6.90
pval = (2*sol - 13)/3
print('Cu average oxidation number', float(2 + pval)); assert round(float(2+pval), 2) == 2.27
print('Cu3+ fraction', float(pval))
# CuO：1 Cu あたり S2O3 1 mol，Cu2O：0
for e in ['2Cu^2+ + 4I- -> 2CuI + I2', 'I2 + 2S2O3^2- -> 2I- + S4O6^2-', 'Cu^3+ + 3I- -> CuI + I2']:
    assert check_equation(e) is None, e
print('equations OK')
# CuO 0.0795 g (1.00e-3 mol) -> 1.00e-3 mol S2O3 -> 50.0 mL of 0.0200 M
print('CuO V mL', float(sp.Rational('0.0795')/(sp.Rational('63.5')+16)/C*1000))
