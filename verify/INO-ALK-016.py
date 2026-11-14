# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem2/out/items/INO-ALK-016/item.json'
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; d = 1
        while d:
            if s[j] == '{': d += 1
            elif s[j] == '}': d -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

eqs = {b for k in ('problem', 'answer', 'solution') for b in ce_blocks(it[k]) if len(re.findall(r'->|<=>', b)) == 1}
for e in sorted(eqs):
    r = check_equation(e); assert r is None, (e, r); print('式OK:', e)

# 図鑑データとの照合
ions = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))
R = ions['results']; cat = {c['id']: c for c in ions['cations']}
assert R['Mg^2+']['NH3']['kind'] == 'ppt' and R['Ca^2+']['NH3']['kind'] == 'none' and R['Ba^2+']['NH3']['kind'] == 'none'
assert R['Mg^2+']['SO4']['kind'] == 'none' and R['Ca^2+']['SO4']['kind'] == 'ppt' and R['Ba^2+']['SO4']['kind'] == 'ppt'
assert cat['Ca^2+']['flame'] == '橙赤' and cat['Ba^2+']['flame'] == '黄緑' and cat['Mg^2+']['flame'] is None
print('図鑑と一致')

# 溶解度積（文献値）で 0.1 mol/L 溶液と試薬を等量混合（各 0.05 mol/L）したときの沈殿の有無
from math import sqrt
Kb = 1.8e-5; c = 0.05
OH = sqrt(Kb * c)                       # アンモニア水 0.1 mol/L 等量混合
Ksp_OH = {'Mg': 5.6e-12, 'Ca': 5.5e-6, 'Ba': 5.0e-3}
Ksp_SO4 = {'Mg': None, 'Ca': 4.9e-5, 'Ba': 1.1e-10}
ppt_NH3 = [m for m in 'Mg Ca Ba'.split() if c * OH**2 > Ksp_OH[m]]
ppt_SO4 = [m for m in 'Mg Ca Ba'.split() if Ksp_SO4[m] and c * c > Ksp_SO4[m]]
print('NH3 で沈殿:', ppt_NH3, ' 希硫酸で沈殿:', ppt_SO4)
assert ppt_NH3 == ['Mg'] and ppt_SO4 == ['Ca', 'Ba']
print('すべて OK')
