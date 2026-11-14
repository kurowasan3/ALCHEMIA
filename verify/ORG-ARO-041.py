# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem2/out/items/ORG-ARO-041/item.json'
it = json.load(open(ITEM, encoding='utf-8'))


def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0:
            return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i + 4:j - 1]); i = j


eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
# 小問数と答えの番号の一致
nq = it['problem'].count('\\item')
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), (nq, na)
print('小問', nq, '答え', na)


def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
from rdkit import Chem
from itertools import combinations
AM = {'H': 1.0, 'C': 12, 'O': 16, 'Br': 80}


def mw(m):
    m = Chem.AddHs(m)
    return sum(AM[a.GetSymbol()] for a in m.GetAtoms())


def free_op(smi):
    m = Chem.MolFromSmiles(smi)
    c1 = [a for a in m.GetAtoms() if a.GetSymbol() == 'C' and a.GetIsAromatic() and any(n.GetSymbol() == 'O' for n in a.GetNeighbors())][0]
    ring = [r for r in m.GetRingInfo().AtomRings() if c1.GetIdx() in r][0]
    k = ring.index(c1.GetIdx())
    pos = [ring[(k + d) % 6] for d in (1, 3, 5)]   # オルト・パラ・オルト
    return m, [i for i in pos if m.GetAtomWithIdx(i).GetTotalNumHs() == 1]


def brominate(m, idxs):
    rw = Chem.RWMol(m)
    for i in idxs:
        j = rw.AddAtom(Chem.Atom(35)); rw.AddBond(i, j, Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(i).SetNoImplicit(False)
    mm = rw.GetMol(); Chem.SanitizeMol(mm)
    return mm


cres = {'phenol': 'Oc1ccccc1', 'o': 'Oc1ccccc1C', 'm': 'Oc1cccc(C)c1', 'p': 'Oc1ccc(C)cc1'}
nbr = {}; mono = {}; prod = {}
for k, s in cres.items():
    m, f = free_op(s)
    nbr[k] = len(f)
    prod[k] = mw(brominate(m, f))
    mono[k] = len({Chem.MolToSmiles(brominate(m, [i])) for i in f})
print('Br数', nbr, '生成物式量', prod, '一臭素置換体', mono)
assert nbr == {'phenol': 3, 'o': 2, 'm': 3, 'p': 2}
assert mono['o'] == 2 and mono['p'] == 1
assert mw(Chem.MolFromSmiles(cres['o'])) == 108
n = 2.16 / 108
assert abs(n - 0.0200) < 1e-9
assert abs(prod['m'] * n - 6.90) < 1e-6 and abs(prod['p'] * n - 5.32) < 1e-6 and abs(prod['o'] * n - 5.32) < 1e-6
# 質量から Br の数を逆算
for mass, want in ((6.90, 3), (5.32, 2)):
    assert round((mass / n - 108) / 79) == want and abs((mass / n - 108) / 79 - want) < 1e-6
nph = 0.940 / mw(Chem.MolFromSmiles('Oc1ccccc1'))
assert sig(3 * nph, 3) == 0.0300
print('ORG-ARO-041 検証OK')
