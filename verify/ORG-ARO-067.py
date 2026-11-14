# ORG-ARO-067 検証：テトラメチルベンゼン → 二無水物 → エタノールによる開環の異性体数
# --- 共通：item.json の \ce{} 反応式のつり合い検査（単独で動く） ---
import json, re, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
def _ce_blocks(s):
    out, i = [], 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0: return out
        k, d = j + 4, 1
        while d:
            if s[k] == '{': d += 1
            elif s[k] == '}': d -= 1
            k += 1
        out.append(s[j+4:k-1]); i = k
def check_item_equations(path):
    it = json.load(open(path, encoding='utf-8'))
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for ce in _ce_blocks(it[key]):
            arrows = re.findall(r'<=>|->|<-', ce)
            if len(arrows) != 1 or '\\boxed' in ce or '\\text' in ce or '\\mathrm' in ce:
                continue
            r = check_equation(ce)
            assert r is None, (key, ce, r)
            n += 1
    print('反応式', n, '本 OK')
    return it
from itertools import combinations, product
from rdkit import Chem
check_item_equations('/home/claude/chemx/todai/out/items/ORG-ARO-067/item.json')
def build(subs, anh=None):
    m = Chem.RWMol(Chem.MolFromSmiles('c1ccccc1'))
    for p, g in subs.items():
        fr = Chem.MolFromSmiles(g); i = m.GetNumAtoms(); m.InsertMol(fr); m.AddBond(p-1, i, Chem.BondType.SINGLE)
    if anh:
        p, q = anh
        c1 = m.AddAtom(Chem.Atom(6)); o1 = m.AddAtom(Chem.Atom(8)); c2 = m.AddAtom(Chem.Atom(6)); oa = m.AddAtom(Chem.Atom(8)); ob = m.AddAtom(Chem.Atom(8))
        m.AddBond(p-1, c1, Chem.BondType.SINGLE); m.AddBond(c1, oa, Chem.BondType.DOUBLE); m.AddBond(c1, o1, Chem.BondType.SINGLE)
        m.AddBond(o1, c2, Chem.BondType.SINGLE); m.AddBond(c2, ob, Chem.BondType.DOUBLE); m.AddBond(c2, q-1, Chem.BondType.SINGLE)
    mol = m.GetMol(); Chem.SanitizeMol(mol); return Chem.MolToSmiles(mol)
tmb = {build({p: 'C' for p in c}) for c in combinations(range(1, 7), 4)}
print('(1)', len(tmb)); assert len(tmb) == 3
adj = lambda a, b: (b - a) % 6 in (1, 5)
maxanh = {}
for pos in [(1, 2, 3, 4), (1, 2, 3, 5), (1, 2, 4, 5)]:
    pairs = [pr for pr in combinations(pos, 2) if adj(*pr)]
    best = max(k for k in range(0, 3) for sel in combinations(pairs, k) if len({x for pr in sel for x in pr}) == 2*k)
    maxanh[pos] = best
print('(2)', maxanh); assert list(maxanh.values()) == [2, 1, 2]
def both_open(r1, r2):
    S = set()
    for c1, c2 in product((0, 1), repeat=2):
        S.add(build({r1[c1]: 'C(=O)OCC', r1[1-c1]: 'C(=O)O', r2[c2]: 'C(=O)OCC', r2[1-c2]: 'C(=O)O'}))
    return S
def one_open(r1, r2):
    S = set()
    for op, kp in [(r1, r2), (r2, r1)]:
        for c in (0, 1):
            S.add(build({op[c]: 'C(=O)OCC', op[1-c]: 'C(=O)O'}, kp))
    return S
n1245 = len(both_open((1, 2), (4, 5))); n1234 = len(both_open((1, 2), (3, 4)))
m1245 = len(one_open((1, 2), (4, 5))); m1234 = len(one_open((1, 2), (3, 4)))
print('(3)', n1245, n1234, '(4)', m1245, m1234)
assert (n1245, n1234, m1245, m1234) == (2, 3, 1, 2)
print('OK')
