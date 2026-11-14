ITEM = '/home/claude/chem2/out/items/ORG-ALC-012/item.json'
# 検証スクリプト（グループ g08）  python3 -I この.py
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

def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return eqs

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def check_parts():
    n_items = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, n_items + 1)], (n_items, nums)
    print('小問の数OK:', n_items)


check_parts()
import itertools
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
RDLogger.DisableLog('rdApp.*')
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
found = set(); atoms = [6, 6, 6, 8]; pairs = list(itertools.combinations(range(4), 2))
for orders in itertools.product(range(2), repeat=6):
    rw = Chem.RWMol()
    for z in atoms: rw.AddAtom(Chem.Atom(z))
    for (a, b), o in zip(pairs, orders):
        if o: rw.AddBond(a, b, Chem.BondType.SINGLE)
    m = rw.GetMol()
    try: Chem.SanitizeMol(m)
    except Exception: continue
    if len(Chem.GetMolFrags(m)) == 1 and not m.GetRingInfo().NumRings() and CalcMolFormula(m) == 'C3H8O': found.add(Chem.MolToSmiles(m))
assert found == {can('CCCO'), can('CC(C)O'), can('CCOC')}, found
bp = {'CCCO': 97, 'CC(C)O': 82, 'CCOC': 7, 'CCC=O': 48, 'CC(C)=O': 56, 'CCC(=O)O': 141}
ox1 = {'CCCO': 'CCC=O', 'CC(C)O': 'CC(C)=O', 'CCOC': 'CCOC'}
ox2 = {'CCC=O': 'CCC(=O)O', 'CC(C)=O': 'CC(C)=O', 'CCOC': 'CCOC'}
X = sorted(ox1, key=lambda s: -bp[s]); assert X == ['CCCO', 'CC(C)O', 'CCOC']
Y = [ox1[x] for x in X]; Z = [ox2[y] for y in Y]
assert all(bp[y] < bp[x] for x, y in zip(X, Y) if x != y)
assert [y != z for y, z in zip(Y, Z)] == [True, False, False] and bp[Z[0]] > bp[X[0]]
TX, TY, TZ = max(bp[s] for s in X), max(bp[s] for s in Y), max(bp[s] for s in Z)
assert TZ > TX > TY
assert len(set(X + Y + Z)) == 6

print('全検証OK:', ITEM)
