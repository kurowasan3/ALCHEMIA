# THE-ELC-043 検証：電池の理論エネルギー密度の比較
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
F = 9.65e4
H, Li, O, Mn, Ni, Zn, Ag, Cd = 1.0, 6.9, 16, 54.9, 58.7, 65.4, 108, 112
cells = {
 'NiCd': (2, 1.30, Cd + 2*(Ni+2*O+H) + 2*(2*H+O)),
 'AgO':  (2, 1.55, Zn + 2*Ag + O),
 'LiMn': (1, 3.00, Li + Mn + 2*O),
}
res = {}
for k, (n, E, m) in cells.items():
    w = n*F*E/m
    res[k] = w
    print(k, 'mass', m, 'J', n*F*E, 'J/g', w, 'g per mol e-', m/n)
assert round(res['NiCd']) == 757
assert round(res['AgO'], -1) == 1010
assert round(res['LiMn'], -1) == 3090
# 水を除いたとき
print('NiCd without water', 2*F*1.30/(Cd + 2*(Ni+2*O+H)))
for eq in ['Cd + 2NiO(OH) + 2H2O -> Cd(OH)2 + 2Ni(OH)2', 'Zn + Ag2O -> ZnO + 2Ag',
           'Li + MnO2 -> LiMnO2', '2Li + 2H2O -> 2LiOH + H2',
           'Cd + 2OH- -> Cd(OH)2 + 2e-', 'NiO(OH) + H2O + e- -> Ni(OH)2 + OH-']:
    r = check_equation(eq); print(eq, r); assert r is None
print('OK')
