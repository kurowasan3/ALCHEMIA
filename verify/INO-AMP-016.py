# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q
ITEM = '/home/claude/chem2/out/items/INO-AMP-016/item.json'
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

Hf = {'Al2O3': -1676, 'Fe2O3': -824, 'Cr2O3': -1140, 'Co3O4': -891, 'MgO': -602, 'CaO': -635}
nO = {'Al2O3': 3, 'Fe2O3': 3, 'Cr2O3': 3, 'Co3O4': 4, 'MgO': 1, 'CaO': 1}
perO = {k: Q(v, nO[k]) for k, v in Hf.items()}
for k, v in perO.items(): print(k, round(float(v), 1))
ok = [k for k in ['Fe2O3', 'Cr2O3', 'Co3O4', 'MgO', 'CaO'] if perO['Al2O3'] - perO[k] < 0]
print('Al で還元できる:', ok); assert ok == ['Fe2O3', 'Cr2O3', 'Co3O4']
dFe = Hf['Al2O3'] - Hf['Fe2O3']; dCo = 4*Hf['Al2O3'] - 3*Hf['Co3O4']; dMg = Hf['Al2O3'] - 3*Hf['MgO']; dCr = Hf['Al2O3'] - Hf['Cr2O3']
print(dFe, dCo, dMg, dCr); assert (dFe, dCo, dMg, dCr) == (-852, -4031, 130, -536)
# O 1 mol あたりの比較と全体の値の一致（Co3O4：O 12 mol）
assert 12*(perO['Al2O3'] - perO['Co3O4']) == dCo
qFe = Q(852, 160 + 2*27); qCr = Q(536, 152 + 2*27)
print(round(float(qFe), 3), round(float(qCr), 3)); assert round(float(qFe), 2) == 3.98 and round(float(qCr), 2) == 2.60
print('すべて OK')
