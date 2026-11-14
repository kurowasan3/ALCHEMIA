# INO-PPT-025 検証：亜鉛による金属の析出と相互の混合による沈殿，亜鉛板の質量変化
import sys, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, nsimplify
D = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))['results']
# イオン化傾向（高校の系列）
series = ['Li', 'K', 'Ca', 'Na', 'Mg', 'Al', 'Zn', 'Fe', 'Ni', 'Sn', 'Pb', 'H', 'Cu', 'Hg', 'Ag', 'Pt', 'Au']
def deposits(metal_ion):
    return series.index(metal_ion) > series.index('Zn')
cand = {'AgNO3': 'Ag', 'Cu(NO3)2': 'Cu', 'Pb(NO3)2': 'Pb', 'NaCl': 'Na', 'MgSO4': 'Mg'}
dep = {k: deposits(v) for k, v in cand.items()}; print(dep)
assert [k for k, v in dep.items() if v] == ['AgNO3', 'Cu(NO3)2', 'Pb(NO3)2']
# 混合で生じる沈殿（図鑑）
assert D['Ag+']['Cl']['f'] == 'AgCl' and D['Ag+']['NH3-ex']['kind'] == 'dissolve'   # B + D → AgCl, NH3 に溶ける
assert D['Pb^2+']['Cl']['f'] == 'PbCl2' and '熱水' in D['Pb^2+']['Cl']['note']       # C + D → PbCl2, 熱水に溶ける
assert D['Pb^2+']['SO4']['f'] == 'PbSO4'                                            # C + E → PbSO4
assert D['Cu^2+']['Cl']['kind'] == 'none' and D['Cu^2+']['SO4']['kind'] == 'none'    # A は D, E と沈殿しない
assert D['Mg^2+']['Cl']['kind'] == 'none'
# 質量変化
Zn, Cu, Ag = R('65.4'), R('63.5'), R(108)
n = R(10, 100) * R(100, 1000)   # 0.010 mol
dA = n * Cu - n * Zn            # Zn + Cu2+ -> Zn2+ + Cu
dB = n * Ag - n / 2 * Zn        # Zn + 2Ag+ -> Zn2+ + 2Ag
print('A:', float(dA), 'B:', float(dB))
assert round(float(dA), 3) == -0.019 and round(float(dB), 2) == 0.75
for eq in ['Zn + Cu^2+ -> Zn^2+ + Cu', 'Zn + 2Ag+ -> Zn^2+ + 2Ag', 'Zn + Pb^2+ -> Zn^2+ + Pb', 'Ag+ + Cl- -> AgCl', 'AgCl + 2NH3 -> [Ag(NH3)2]+ + Cl-']:
    assert check_equation(eq) is None, eq
print('INO-PPT-025 OK')
