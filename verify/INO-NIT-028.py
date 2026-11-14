# INO-NIT-028 検証：リン酸カルシウムと塩酸の量的関係
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
eqs = ['Ca3(PO4)2 + 2HCl -> 2CaHPO4 + CaCl2',
       'Ca3(PO4)2 + 4HCl -> Ca(H2PO4)2 + 2CaCl2',
       'Ca3(PO4)2 + 6HCl -> 2H3PO4 + 3CaCl2',
       '2CaHPO4 + 2HCl -> Ca(H2PO4)2 + CaCl2',
       'PO4^3- + H+ -> HPO4^2-']
for e in eqs:
    r = check_equation(e); print(e, r); assert r is None
# 0.0100 mol Ca3(PO4)2 + 0.0300 mol HCl
from fractions import Fraction as F
n0 = F(1, 100); h = F(3, 100)
# 第1段：2 当量で CaHPO4 2n0 + CaCl2 n0
CaHPO4 = 2*n0; CaCl2 = n0; h -= 2*n0
# 第2段：2CaHPO4 + 2HCl -> Ca(H2PO4)2 + CaCl2
x = h  # HCl 消費
CaHPO4 -= x; CaH2 = x/2; CaCl2 += x/2
print(float(CaHPO4), float(CaH2), float(CaCl2))
assert (CaHPO4, CaH2, CaCl2) == (F(1,100), F(1,200), F(3,200))
# 原子の保存
assert CaHPO4 + CaH2 + CaCl2 == 3*n0          # Ca
assert CaHPO4 + 2*CaH2 == 2*n0                # P
assert 2*CaCl2 == F(3,100)                    # Cl
# 全部溶かす（H3PO4 まで）に要する HCl：6 n0 = 0.0600 mol → 2.0 mol/L なら 30 mL
print('HCl for complete =', float(6*n0), 'mol ->', float(6*n0/2*1000), 'mL of 2.0 mol/L')
print('OK')
