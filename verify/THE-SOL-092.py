# THE-SOL-092 検証：直線で近似した蒸気圧曲線と蒸気圧降下から沸点上昇・モル沸点上昇を導く
from sympy import Rational as Q, symbols, solve, nsimplify
P0 = 1.013e5          # 大気圧 Pa
A_bz = 3.10e3         # ベンゼンの沸点付近の dp/dT (Pa/K)
M_bz = 0.078          # kg/mol
m = 0.200             # mol/kg
x2 = m * M_bz
x2_exact = m / (m + 1 / M_bz)
dp = x2 * P0
dT = dp / A_bz
Kb = M_bz * P0 / A_bz
print('x2', x2, 'exact', x2_exact, 'dp', dp, 'dT', dT, 'Kb', Kb, 'Kb*m', Kb * m)
assert f'{x2:.3g}' == '0.0156' and f'{dp:.3g}' == '1.58e+03' and f'{dT:.2g}' == '0.51' and f'{Kb:.2g}' == '2.5'
# 溶液の沸点：P_sol(T) = (1-x2)(P0 + A (T - Tb)) = P0 を厳密に解いた場合との差
T = symbols('T')
Tb = 353.2
sol = solve((1 - x2) * (P0 + A_bz * (T - Tb)) - P0, T)[0]
print('exact dT (Raoult on line)', float(sol - Tb))
assert abs(float(sol - Tb) - dT) / dT < 0.02
# 水
A_w = 3.62e3; M_w = 0.018
Kb_w = M_w * P0 / A_w
print('Kb water', Kb_w)
assert f'{Kb_w:.2g}' == '0.5'
print('ratio Kb', Kb / Kb_w, 'M ratio', M_bz / M_w, 'slope ratio', A_w / A_bz)
print('OK')
