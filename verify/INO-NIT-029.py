# INO-NIT-029 検証：二リン酸ナトリウムによる硬水の軟化とリンの排水基準
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['2Na2HPO4 -> Na4P2O7 + H2O','2Ca^2+ + P2O7^4- -> Ca2P2O7','2Mg^2+ + P2O7^4- -> Mg2P2O7']:
    print(e, check_equation(e) or 'OK')
H,O,Na,Mg,P,Ca,C = 1.0,16,23,24.3,31,40,12
M_Na4P2O7 = 4*Na+2*P+7*O; M_CaCO3 = Ca+C+3*O
print('M Na4P2O7', M_Na4P2O7, 'CaCO3', M_CaCO3)
# 水の組成 Ca 32.0 mg/L, Mg 9.72 mg/L
nCa = 32.0e-3/Ca; nMg = 9.72e-3/Mg
hard = (nCa+nMg)*M_CaCO3*1000
print('mmol/L Ca,Mg', nCa*1e3, nMg*1e3, 'hardness mg/L', hard)
V = 50.0
# P 上限 16 mg/L
nP = 16e-3*V/P; mass = nP/2*M_Na4P2O7
print('nP', nP, 'mass Na4P2O7 g', mass)
nM = (nCa+nMg)*V; nP2O7 = nP/2
removed = 2*nP2O7
rest = nM - removed
print('M2+ total mmol', nM*1e3, 'removed', removed*1e3, 'rest', rest*1e3, 'residual hardness mg/L', rest/V*M_CaCO3*1000)
print('removed fraction', removed/nM)
# 全部除くのに必要な P 濃度
needP = (nCa+nMg)*1*P*1000   # mol P per L = nM per L (P2O7 binds 2 M, has 2 P)
print('P needed mg/L', needP)
# 基準内で完全に軟化できる最大硬度
print('max hardness (CaCO3 mg/L) at P 16 mg/L', 16e-3/P*M_CaCO3*1000)
