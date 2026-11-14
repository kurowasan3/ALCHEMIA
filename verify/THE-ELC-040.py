"""THE-ELC-040 検証：2室の電気分解での電極反応の切り替わりと pH（SymPy で独立計算）"""
import sys, json, re
from sympy import Rational as R, log, N
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

nOH0 = R(500, 10000)   # 陽極室 NaOH 0.0500 mol（1.00 L）
nCu0 = R(400, 10000)   # 陰極室 CuSO4 0.0400 mol
Kw = R(1, 10**14)
def pH_anode(x):
    if x < nOH0: return 14 + log(nOH0 - x, 10)   # [OH-] = 0.0500 - x
    return -log(x - nOH0, 10)                      # [H+] = x - 0.0500
def pH_cath(x):
    xs = 2 * nCu0                                  # Cu2+ 還元に要する電子 0.0800 mol
    return 7 if x <= xs else 14 + log(x - xs, 10)
vals = {x: round(float(N(pH_anode(R(x)))), 2) for x in ['0', '3/100', '49/1000', '51/1000', '7/100']}
print('陽極室 pH', vals)
assert vals == {'0': 12.70, '3/100': 12.30, '49/1000': 11.00, '51/1000': 3.00, '7/100': 1.70}
print('気体発生開始（陰極）', float(2 * nCu0))
assert 2 * nCu0 == R(8, 100)
pc = round(float(N(pH_cath(R(1, 10)))), 2); print('陰極室 pH (x=0.100)', pc); assert pc == 12.30
# 0.100 mol 流したときの気体
O2 = R(1, 10) / 4; H2 = (R(1, 10) - 2 * nCu0) / 2; print('O2', O2, 'H2', H2)
assert O2 == R(25, 1000) and H2 == R(1, 100)
# log2 = 0.30 を使った値との一致（12.70, 12.30, 1.70）
assert abs(float(14 + log(R(5,100),10)) - 12.70) < 0.005
# 反応式のつり合い
for eq in ['4OH- -> O2 + 2H2O + 4e-', '2H2O -> O2 + 4H+ + 4e-', 'Cu^2+ + 2e- -> Cu', '2H2O + 2e- -> H2 + 2OH-']:
    r = check_equation(eq); print(eq, r or 'OK'); assert r is None
print('THE-ELC-040 all OK')
