# INO-PPT-028 検証：液性・共存物で沈殿の有無が変わる例
import sys, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
R = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))['results']
# (1)(2)(3) CO2 と Ca2+：弱酸 H2CO3 は強酸 HCl を遊離できない（(1) 沈殿なし），塩基性なら CaCO3，過剰 CO2 で Ca(HCO3)2
assert R['Ca^2+']['CO3']['f'] == 'CaCO3' and R['Ca^2+']['CO3']['color'] == '白'
# 平衡の目安：CaCO3 + 2H+ -> Ca2+ + H2CO3 の K = Ksp/(Ka1 Ka2) が大きい -> 酸性（中性の CaCl2 + CO2 飽和）では沈殿しない
Ksp, Ka1, Ka2 = 4.8e-9, 4.5e-7, 4.7e-11
K = Ksp / (Ka1 * Ka2); print('K(CaCO3 dissolve in acid)=%.1e' % K); assert K > 1e7
# CO2 飽和水 (H2CO3 ~0.033 M, pH ~ 3.9) で [CO3 2-] ~ Ka2 = 4.7e-11 -> [Ca2+]=0.1 との積 4.7e-12 << Ksp
assert 0.1 * Ka2 < Ksp
# (4) Fe3+ + H2S（酸性）-> S（白濁），FeS は沈殿しない
assert R['Fe^3+']['H2S-a']['f'] == 'S'
# (5) Zn2+：酸性 H2S で沈殿なし，塩基性で ZnS（白）
assert R['Zn^2+']['H2S-a']['kind'] == 'none' and R['Zn^2+']['H2S-b']['f'] == 'ZnS' and R['Zn^2+']['H2S-b']['color'] == '白'
# Zn2+ に NH3 を過剰に加えると [Zn(NH3)4]2+ だが，H2S を通じると ZnS が沈殿（ZnS は非常に溶けにくい）
eqs = ['CaCl2 + 2NH3 + CO2 + H2O -> CaCO3 + 2NH4Cl', 'CaCO3 + CO2 + H2O -> Ca(HCO3)2',
       '2Fe^3+ + H2S -> 2Fe^2+ + S + 2H+', 'Zn^2+ + S^2- -> ZnS', '[Zn(NH3)4]^2+ + H2S -> ZnS + 2NH4+ + 2NH3']
for e in eqs:
    assert check_equation(e) is None, e
print('INO-PPT-028 OK')
