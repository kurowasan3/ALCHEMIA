# THE-ACB-031 検証：硫酸アンモニウムと塩化アンモニウムの混合物（NH3 逆滴定＋BaSO4 重量分析）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from sympy import symbols, solve, Rational as R
from chemlib import check_equation
H, N, O, S, Cl, Ba = 1.0, 14, 16, 32, 35.5, 137
M_AS = 2*(N+4*H) + S + 4*O      # (NH4)2SO4
M_AC = N + 4*H + Cl             # NH4Cl
M_B = Ba + S + 4*O              # BaSO4
print('M', M_AS, M_AC, M_B)
# 分取 10.0 mL（全体 100 mL の 1/10）
H_acid = 2*0.100*25.0e-3        # 吸収用 H2SO4 の H+
OH_back = 0.100*18.0e-3         # 逆滴定の NaOH
n_NH3 = H_acid - OH_back
print('NH3 in aliquot', n_NH3)
n_SO4 = 0.233 / M_B
print('SO4 in aliquot', n_SO4)
x, y = symbols('x y')
sol = solve([2*x + y - 10*n_NH3, x - 10*n_SO4], [x, y])
print(sol, 'mass AS', float(sol[x]*M_AS), 'mass AC', float(sol[y]*M_AC))
assert abs(float(sol[x]) - 0.0100) < 1e-6 and abs(float(sol[y]) - 0.0120) < 1e-6
assert round(float(sol[x]*M_AS), 2) == 1.32 and round(float(sol[y]*M_AC), 3) == 0.642
tot = float(sol[x]*M_AS + sol[y]*M_AC)
print('total', tot, 'NH4Cl %', float(sol[y]*M_AC)/tot*100)
# 吸収液の硫酸から生じうる BaSO4
print('BaSO4 from absorber', 0.100*25.0e-3*M_B)
# 吸収液に NH3 が全部吸収できる（H+ が過剰）こと
assert H_acid > n_NH3
for eq in ['NH4+ + OH- -> NH3 + H2O', '(NH4)2SO4 + 2NaOH -> Na2SO4 + 2NH3 + 2H2O',
           'NH4Cl + NaOH -> NaCl + NH3 + H2O', '2NH3 + H2SO4 -> (NH4)2SO4',
           'H2SO4 + 2NaOH -> Na2SO4 + 2H2O', 'Ba^2+ + SO4^2- -> BaSO4']:
    r = check_equation(eq); print(eq, r); assert r is None
print('OK')
