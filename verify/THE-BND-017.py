# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-017/item.json'
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
lp = {'N': 1, 'O': 2, 'F': 3}
def hb(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    h = sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'H' and a.GetNeighbors()[0].GetSymbol() in lp)
    l = sum(lp.get(a.GetSymbol(), 0) for a in m.GetAtoms() if a.GetSymbol() in lp)
    return h, l, min(h, l)
r = {k: hb(s) for k, s in {'H2O': 'O', 'HF': 'F', 'NH3': 'N', 'DME': 'COC'}.items()}
print(r)
assert r['H2O'][2] == 2 and r['HF'][2] == 1 and r['NH3'][2] == 1 and r['DME'] == (0, 2, 0)
bp = {'H2O': 100, 'HF': 20, 'NH3': -33}
assert bp['H2O'] > bp['HF'] > bp['NH3']

print('全検証OK:', ITEM)
