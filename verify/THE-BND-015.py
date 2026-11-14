# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-015/item.json'
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and not re.search(r'_n|\bn[A-Z]|R', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
mw = {k: MolWt(Chem.MolFromSmiles(s)) for k, s in {'F2': 'FF', 'Cl2': 'ClCl', 'N2': 'N#N', 'CO': '[C-]#[O+]', 'DME': 'COC', 'EtOH': 'CCO'}.items()}
print(mw)
assert abs(mw['N2'] - mw['CO']) < 0.01 and abs(mw['DME'] - mw['EtOH']) < 1e-6
bp = {'F2': -188, 'Cl2': -34, 'N2': -196, 'CO': -192, 'DME': -24, 'EtOH': 78}
mp = {'diamond': 3550, 'Si': 1410, 'NaF': 993, 'MgO': 2826}
ans = ['Cl2', 'CO', 'EtOH', 'diamond', 'MgO']
pairs = [('F2', 'Cl2', bp), ('N2', 'CO', bp), ('DME', 'EtOH', bp), ('diamond', 'Si', mp), ('NaF', 'MgO', mp)]
for (a, b, t), w in zip(pairs, ans):
    hi = a if t[a] > t[b] else b
    print(a, b, '→', hi); assert hi == w
r = {'Na+': 0.102, 'F-': 0.133, 'Mg2+': 0.072, 'O2-': 0.140}
print('NaF 距離', r['Na+'] + r['F-'], ' MgO 距離', r['Mg2+'] + r['O2-'])

print('全検証OK:', ITEM)
