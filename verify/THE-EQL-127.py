# THE-EQL-127 検証：速度定数の比による平衡定数と窒素原子の保存（N2O5 の分解 → NO2/N2O4 平衡）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from sympy import Rational as Q, symbols, solve, nsimplify
from chemlib import check_equation
for eq in ['2N2O5 -> 4NO2 + O2', '2NO2 <=> N2O4']:
    assert check_equation(eq) is None, eq
k1 = 4.0e8; k2 = 2.5e6
K = k1 / k2
assert abs(K - 160) < 1e-9
r = 2.50; V = 4.00
no2 = 1 / (r * K); n2o4 = no2 / r
assert abs(n2o4 - K * no2**2) < 1e-15
N = V * (no2 + 2 * n2o4)            # 窒素原子 mol
n2o5 = N / 2
o2 = n2o5 / 2
print('K', K, '[NO2]', no2, '[N2O4]', n2o4, 'N2O5', n2o5, 'O2', o2)
assert abs(no2 - 2.5e-3) < 1e-12 and abs(n2o4 - 1.0e-3) < 1e-12 and abs(n2o5 - 9.0e-3) < 1e-12
R = 8.31e3; T = 300
P = (V * (no2 + n2o4) + o2) * R * T / V
print('P', P)
assert f'{P:.3g}' == '1.15e+04'
# 2.00 L に圧縮したときの比
x = symbols('x', positive=True)
Nc = N / 2.00
sol = [s for s in solve(2 * K * x**2 + x - Nc, x) if s > 0][0]
ratio = float(1 / (K * sol))
print('ratio after compression', ratio)
assert ratio < 2.50 and f'{ratio:.3g}' == '1.58'
print('OK')
