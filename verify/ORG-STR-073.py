# 検証スクリプト（グループ g03）  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/ORG-STR-073/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from collections import defaultdict

def nmr(smi):
    """等価な H の組（対称性クラス）ごとに (H数, 隣接炭素上の H 数 n → n+1 本) を返す"""
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    rank = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    groups = defaultdict(list)
    for a in m.GetAtoms():
        if a.GetSymbol() == 'H':
            groups[rank[a.GetIdx()]].append(a)
    out = []
    for g in groups.values():
        c = g[0].GetNeighbors()[0]
        nb = 0
        for x in c.GetNeighbors():
            if x.GetSymbol() == 'C':
                nb += sum(1 for y in x.GetNeighbors() if y.GetSymbol() == 'H')
        env = 'O' if any(x.GetSymbol() == 'O' for x in c.GetNeighbors()) else (
              'CO' if any(x.GetSymbol() == 'C' and any(b.GetBondTypeAsDouble() == 2 for b in x.GetBonds()) for x in c.GetNeighbors()) else 'alkyl')
        out.append((len(g), nb + 1, env))
    return sorted(out)

cases = {
 'P': ('CCC(=O)OCC', [(3, 3, 'alkyl'), (3, 3, 'alkyl'), (2, 4, 'CO'), (2, 4, 'O')]),
 'Q': ('CC(=O)OC(C)C', [(6, 2, 'alkyl'), (3, 1, 'CO'), (1, 7, 'O')]),
 'R': ('CC(C)C(=O)OC', [(6, 2, 'alkyl'), (1, 7, 'CO'), (3, 1, 'O')]),
 'propyl acetate': ('CC(=O)OCCC', [(3, 3, 'alkyl'), (2, 6, 'alkyl'), (3, 1, 'CO'), (2, 3, 'O')]),
}
for k, (smi, exp) in cases.items():
    m = Chem.MolFromSmiles(smi)
    f = CalcMolFormula(m)
    got = nmr(smi)
    print(k, f, got)
    assert f == 'C5H10O2'
    assert got == sorted(exp), (k, got)
# 化学シフトが与えた範囲に入るか（O-CH: 3.5-5.2, C(=O)-CH: 2.0-2.7, alkyl 0.8-1.8）
rng = {'O': (3.5, 5.2), 'CO': (2.0, 2.7), 'alkyl': (0.8, 1.8)}
data = {'P': [(1.14, 'alkyl'), (1.26, 'alkyl'), (2.32, 'CO'), (4.12, 'O')],
        'Q': [(1.23, 'alkyl'), (2.02, 'CO'), (5.00, 'O')],
        'R': [(1.17, 'alkyl'), (2.55, 'CO'), (3.68, 'O')]}
for k, v in data.items():
    for d, e in v:
        assert rng[e][0] <= d <= rng[e][1], (k, d, e)
# Q と R は分裂・強度だけでは同じ
strip = lambda L: sorted((a, b) for a, b, _ in L)
assert strip(nmr('CC(=O)OC(C)C')) == strip(nmr('CC(C)C(=O)OC'))
print('OK')
