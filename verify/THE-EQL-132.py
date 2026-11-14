# THE-EQL-132 検証：N2O3 の分解・NO の酸化（O2 が過不足どちらか）・NO2 の二量化平衡の連成
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from sympy import symbols, sqrt, nsolve, Rational as Q
from chemlib import check_equation
for eq in ['N2O3 -> NO + NO2', '2NO + O2 -> 2NO2', '2NO2 <=> N2O4']:
    assert check_equation(eq) is None, eq
K = 100.0; V = 1.00; R = 8.31e3; T = 300.0; b = 5.00e-3
def nsys(N):   # N: NO2/N2O4 系の窒素原子 mol（V=1 L なので濃度と同じ）
    x = (-1 + (1 + 8 * K * N / V) ** 0.5) / (4 * K / V)   # 2K x^2 + x - N = 0（V=1）
    y = K * x * x / V
    return x, y
def ntot(a):
    if a <= 2 * b:          # O2 十分：NO はすべて NO2 に
        x, y = nsys(2 * a); return x + y + (b - a / 2)
    else:                   # O2 不足：NO が残る
        x, y = nsys(a + 2 * b); return x + y + (a - 2 * b)
# 単調性
import numpy as np
grid = np.linspace(1e-4, 3e-2, 3000)
vals = [ntot(a) for a in grid]
assert all(v2 > v1 for v1, v2 in zip(vals, vals[1:]))
a_true = 8.0e-3
P_true = ntot(a_true) * R * T / V
P_given = float(f'{P_true:.3g}')
print('P true', P_true, 'given', P_given)
Pb = ntot(2 * b) * R * T / V
print('boundary P (a=2b)', Pb)
assert P_given < Pb
# 逆算（二分法）
lo, hi = 1e-5, 2 * b
for _ in range(200):
    mid = (lo + hi) / 2
    if ntot(mid) * R * T / V < P_given: lo = mid
    else: hi = mid
a = (lo + hi) / 2
print('a', a)
assert f'{a:.2g}' == '0.008'
x, y = nsys(2 * a)
pO2 = (b - a / 2) * R * T / V
print('NO2', x, 'N2O4', y, 'O2 left', b - a / 2, 'pO2', pO2)
# 必要最小の O2
print('O2 min for a', a / 2)
# 場合2 を仮定した解が存在しない（a>2b で ntot>P_given）ことを確認
assert ntot(2 * b + 1e-9) * R * T / V > P_given
print('OK')
