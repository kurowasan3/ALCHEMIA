# INO-OXS-014 検証：陰イオン（CO3^2-, SO3^2-, SO4^2-, Cl-）の識別
import sys, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
# 各候補の観察（高校教科書の標準的事実）
# BaCl2 で沈殿：BaCO3, BaSO3, BaSO4（白）／ Cl-, NO3- は沈殿なし
# 希塩酸：BaCO3 -> CO2（無臭），BaSO3 -> SO2（刺激臭，KMnO4 脱色），BaSO4 は溶けない（強酸の塩）
cand = {
 'Na2CO3': dict(ba=True, hcl='gas', odor=False, kmno4=False, ag=None),
 'Na2SO3': dict(ba=True, hcl='gas', odor=True, kmno4=True, ag=None),
 'Na2SO4': dict(ba=True, hcl='insoluble', odor=None, kmno4=None, ag=None),
 'NaCl':   dict(ba=False, ag='white, dissolves in NH3'),
 'NaNO3':  dict(ba=False, ag=None),
}
obs = {'A': dict(ba=True, hcl='gas', odor=False), 'B': dict(ba=True, hcl='gas', odor=True, kmno4=True),
       'C': dict(ba=True, hcl='insoluble'), 'D': dict(ba=False, ag='white, dissolves in NH3')}
for L, o in obs.items():
    m = [k for k, v in cand.items() if all(v.get(kk) == vv for kk, vv in o.items())]
    print(L, m); assert len(m) == 1
# 石灰水：CO2 も SO2 も白濁（CaCO3, CaSO3 はともに水に溶けにくい） -> 区別できない
eqs = ['BaCl2 + Na2CO3 -> BaCO3 + 2NaCl', 'BaSO3 + 2HCl -> BaCl2 + SO2 + H2O', 'BaCO3 + 2HCl -> BaCl2 + CO2 + H2O',
       'Ca(OH)2 + CO2 -> CaCO3 + H2O', 'Ca(OH)2 + SO2 -> CaSO3 + H2O',
       'SO3^2- + I2 + H2O -> SO4^2- + 2I- + 2H+', 'Ba^2+ + SO4^2- -> BaSO4',
       '5SO2 + 2MnO4- + 2H2O -> 5SO4^2- + 2Mn^2+ + 4H+']
for e in eqs:
    assert check_equation(e) is None, e
D = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))['results']
assert D['Ba^2+']['SO4']['f'] == 'BaSO4' and D['Ba^2+']['CO3']['f'] == 'BaCO3'
assert D['Ag+']['Cl']['f'] == 'AgCl' and D['Ag+']['NH3-ex']['kind'] == 'dissolve'
print('INO-OXS-014 OK')
