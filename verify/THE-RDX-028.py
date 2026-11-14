# THE-RDX-028 検証：ビタミンC（アスコルビン酸）のヨウ素による直接滴定
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
H, C, O = 1.0, 12, 16
M = 6*C + 8*H + 6*O
print('M(C6H8O6)=', M); assert M == 176
nI2 = 0.0100 * 11.4e-3          # 10.0 mL 中
n10 = nI2                        # 1:1
n_tab = n10 * 100/10
mg = n_tab * M * 1000
print('mg per tablet', mg, round(mg, 0)); assert round(mg) == 201
# KMnO4（電子5）で同じ 10.0 mL を滴定：e- = 2*n10 = 5*nMn
nMn = 2*n10/5
V = nMn/0.00200*1000
print('KMnO4 mL', V); assert abs(V-22.8) < 0.05
for eq in ['C6H8O6 -> C6H6O6 + 2H+ + 2e-', 'C6H8O6 + I2 -> C6H6O6 + 2HI',
           '5C6H8O6 + 2MnO4- + 6H+ -> 5C6H6O6 + 2Mn^2+ + 8H2O',
           '2C6H8O6 + O2 -> 2C6H6O6 + 2H2O', 'I2 + 2S2O3^2- -> 2I- + S4O6^2-']:
    r = check_equation(eq); print(eq, r); assert r is None
print('OK')
