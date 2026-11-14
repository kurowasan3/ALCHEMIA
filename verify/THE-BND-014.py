# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-014/item.json'
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
from rdkit.Chem import AllChem
import numpy as np
VAL = {'H': 1, 'C': 4, 'N': 5, 'O': 6}
def pairs(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    shared = sum(b.GetBondTypeAsDouble() for b in m.GetBonds())
    ve = sum(VAL[a.GetSymbol()] for a in m.GetAtoms()) - Chem.GetFormalCharge(m)
    return int(shared), int((ve - 2 * shared) / 2), m.GetNumAtoms()
cand = {'NH3': 'N', 'CO2': 'O=C=O', 'N2': 'N#N', 'HCN': 'C#N', 'CO': '[C-]#[O+]', 'N2O': '[N-]=[N+]=O'}
for k, s in cand.items(): print(k, pairs(s))
assert pairs('N') == (3, 1, 4)
assert pairs('O=C=O') == (4, 4, 3) and pairs('[N-]=[N+]=O')[:2] == (4, 4)
assert pairs('N#N') == (3, 2, 2) and pairs('[C-]#[O+]')[:2] == (3, 2)
assert pairs('C#N') == (4, 1, 3)

print('全検証OK:', ITEM)
