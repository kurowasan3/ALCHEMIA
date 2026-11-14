# 検証スクリプト（グループ g03）  python3 -I この.py

import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

def check_all_ce(it):
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if should_check(b):
                r = check_equation(clean_ce(b))
                assert r is None, (b, r)
                print('式OK:', b); n += 1
    print('検査した反応式:', n)

from rdkit import Chem
from itertools import product
it = json.load(open('/home/claude/chem5/out/items/ORG-STR-038/item.json', encoding='utf-8'))
check_all_ce(it)
# C3H6O の構造異性体の総当たり
atoms = ['C', 'C', 'C', 'O']; pairs = [(i, j) for i in range(4) for j in range(i+1, 4)]
val = {'C': 4, 'O': 2}; bt = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
iso = set()
for orders in product([0, 1, 2, 3], repeat=6):
    deg = [0]*4
    for (i, j), o in zip(pairs, orders): deg[i] += o; deg[j] += o
    if any(deg[k] > val[atoms[k]] for k in range(4)) or sum(val[atoms[k]] - deg[k] for k in range(4)) != 6: continue
    rw = Chem.RWMol()
    for a in atoms: rw.AddAtom(Chem.Atom(a))
    for (i, j), o in zip(pairs, orders):
        if o: rw.AddBond(i, j, bt[o])
    m = rw.GetMol()
    try: Chem.SanitizeMol(m)
    except Exception: continue
    if len(Chem.GetMolFrags(m)) == 1: iso.add(Chem.MolToSmiles(m))
print('構造異性体', len(iso), sorted(iso)); assert len(iso) == 9
enol = Chem.MolFromSmarts('C=C[OX2H]')
stable = sorted(s for s in iso if not Chem.MolFromSmiles(s).HasSubstructMatch(enol))
assert len(stable) == 7
pat = {k: Chem.MolFromSmarts(v) for k, v in {'OH': '[CX4][OX2H]', 'CC': 'C=C', 'CHO': '[CX3H1](=O)', 'I': '[CH3]C(=O)'}.items()}
pat2 = Chem.MolFromSmarts('[CH3][CH1][OH]')
def tests(s):
    m = Chem.MolFromSmiles(s)
    return (m.HasSubstructMatch(pat['OH']), m.HasSubstructMatch(pat['CC']), m.HasSubstructMatch(pat['CHO']),
            m.HasSubstructMatch(pat['I']) or m.HasSubstructMatch(pat2))
T = {s: tests(s) for s in stable}
for s in stable: print(s, T[s])
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
P = [s for s in stable if T[s][0] and not T[s][1]]; assert P == [can('OC1CC1')]
none = sorted(s for s in stable if not any(T[s])); assert set(none) == {can('C1COC1'), can('CC1CO1')}
def nC(s): return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
assert nC('CC1CO1') == 1 and nC('C1COC1') == 0
Q = [s for s in stable if T[s][1] and T[s][0]]; R = [s for s in stable if T[s][1] and not T[s][0]]
assert Q == [can('C=CCO')] and R == [can('C=COC')]
assert nC('BrCC(Br)CO') == 1 and nC('BrCC(Br)OC') == 1
X = [s for s in stable if T[s][2]]; Y = [s for s in stable if T[s][3]]
assert X == [can('CCC=O')] and Y == [can('CC(C)=O')]
print('すべて一致')
