# INO-PPT-027 検証：系統分析の操作の意味（省いたとき・試薬を替えたときの影響）
import sys, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
D = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))
R = D['results']
# 標準の流れ
assert R['Ag+']['Cl']['kind'] == 'ppt' and all(R[i]['Cl']['kind'] == 'none' for i in ['Cu^2+', 'Fe^3+', 'Zn^2+', 'Ca^2+', 'Na+'])
assert R['Cu^2+']['H2S-a']['f'] == 'CuS'
assert R['Fe^3+']['H2S-a']['f'] == 'S' and 'Fe²⁺' in R['Fe^3+']['H2S-a']['note']   # 操作2で Fe3+ は Fe2+ に還元される
assert R['Zn^2+']['H2S-a']['kind'] == 'none' and R['Ca^2+']['H2S-a']['kind'] == 'none'
assert R['Fe^3+']['NH3-ex']['kind'] == 'remain' and R['Zn^2+']['NH3-ex']['kind'] == 'dissolve'
assert R['Zn^2+']['H2S-b']['f'] == 'ZnS' and R['Ca^2+']['H2S-b']['kind'] == 'none'
assert R['Ca^2+']['CO3']['f'] == 'CaCO3'
# (1) H2S を追い出さないと，塩基性にしたとき ZnS が沈殿 C に混じる
assert R['Zn^2+']['H2S-b']['kind'] == 'ppt'
# (2) Fe2+ のまま：酸性 H2S では沈殿しないが，塩基性 H2S（操作5）では FeS（黒）として沈殿 → D に混じる
assert R['Fe^2+']['H2S-b']['f'] == 'FeS' and R['Fe^2+']['H2S-b']['color'] == '黒'
# (3) NaOH 過剰でも Zn は [Zn(OH)4]2- として溶け Fe(OH)3 は残る → 分離は可能だが Na+ を加えてしまう
assert R['Zn^2+']['NaOH-ex']['kind'] == 'dissolve' and R['Fe^3+']['NaOH-ex']['kind'] == 'remain'
# (4) 希硫酸：Ca2+ が CaSO4 として沈殿 A に混じる
assert R['Ca^2+']['SO4']['f'] == 'CaSO4'
print('systematic (図鑑):', [s['group'] for s in D['systematic']])
for eq in ['2Fe^3+ + H2S -> 2Fe^2+ + S + 2H+', 'Zn^2+ + H2S -> ZnS + 2H+', '3Fe^2+ + NO3- + 4H+ -> 3Fe^3+ + NO + 2H2O',
           'Fe^2+ + H2S -> FeS + 2H+']:
    assert check_equation(eq) is None, eq
print('INO-PPT-027 OK')
