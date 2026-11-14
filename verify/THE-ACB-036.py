# THE-ACB-036 検証：二段滴定の途中で CO2 が逃げたときの滴下量の配分と組成の誤差
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['OH- + H+ -> H2O', 'CO3^2- + H+ -> HCO3-', 'HCO3- + H+ -> H2O + CO2', '2HCO3- -> CO3^2- + H2O + CO2']:
    assert check_equation(e) is None, e
from fractions import Fraction as F
a, b, x, c = F(1), F(2), F(1, 5), F(1, 10)     # mmol, mmol, mmol, mol/L
V1, V2 = (a + b)/c, b/c
print('正常', V1, V2); assert (V1, V2) == (30, 20)
# 第1中和点で x の CO2 が逃げる：HCO3- が 2x 減り CO3^2- が x 生じる → 再び無色まで x，第2段は (b-2x+x)
V1p, V2p = (a + b + x)/c, (b - x)/c
print('CO2 損失', V1p, V2p, '合計', V1p + V2p); assert (V1p, V2p) == (32, 18) and V1p + V2p == V1 + V2
bp = V2p*c; ap = V1p*c - bp
print('見かけの NaOH', ap, ' Na2CO3', bp); assert ap == a + 2*x and bp == b - x
print('OK')
