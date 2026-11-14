# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem2/out/items/ORG-ARO-040/item.json'
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
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import combinations
S = {'A': 'C1=CC=CC=C1', 'B': 'C1=CC2C=CC12', 'C': 'C12C3C1C4C2C34', 'D': 'C1=CC1C1C=C1'}
F = Chem.SanitizeFlags


def canon(m):
    m2 = Chem.Mol(m)
    Chem.SanitizeMol(m2, F.SANITIZE_ALL ^ F.SANITIZE_SETAROMATICITY)
    Chem.RemoveStereochemistry(m2)
    return Chem.MolToSmiles(m2, kekuleSmiles=True)


def subst(smi, comb, arom=False):
    m = Chem.MolFromSmiles(smi, sanitize=False)
    rw = Chem.RWMol(m)
    for i in comb:
        j = rw.AddAtom(Chem.Atom(17)); rw.AddBond(i, j, Chem.BondType.SINGLE)
    mm = rw.GetMol()
    return Chem.CanonSmiles(Chem.MolToSmiles(mm)) if arom else canon(mm)


def count(smi, k, arom=False):
    return {subst(smi, c, arom) for c in combinations(range(6), k)}


res = {k: [len(count(s, n)) for n in (1, 2, 3)] for k, s in S.items()}
res['benzene'] = [len(count('c1ccccc1', n, True)) for n in (1, 2, 3)]
print(res)
assert res == {'A': [1, 4, 4], 'B': [2, 6, 6], 'C': [1, 3, 3], 'D': [2, 5, 5], 'benzene': [1, 3, 3]}
# 三角柱の頂点の番号付け：上 1,2,3 / 下 4,5,6（1-4, 2-5, 3-6）
P = Chem.RWMol()
for _ in range(6): P.AddAtom(Chem.Atom(6))
for a, b in [(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)]: P.AddBond(a, b, Chem.BondType.SINGLE)
ps = Chem.MolToSmiles(P.GetMol())
m = P.GetMol(); Chem.SanitizeMol(m)
assert CalcMolFormula(Chem.AddHs(m)) == 'C6H6'
assert Chem.CanonSmiles(ps) == Chem.CanonSmiles(S['C'])
di = {subst(ps, c) for c in [(0,1),(0,3),(0,4)]}
tri = {subst(ps, c) for c in [(0,1,2),(0,1,3),(0,1,5)]}
assert di == count(ps, 2) and tri == count(ps, 3)
print('三角柱の代表 1,2/1,4/1,5 と 1,2,3/1,2,4/1,2,6 で全種類を尽くす')
mc = Chem.MolFromSmiles(S['C'])
assert sum(b.GetBondTypeAsDouble() == 2 for b in mc.GetBonds()) == 0
assert CalcMolFormula(Chem.MolFromSmiles('C1CCCCC1')) == 'C6H12'
print('ORG-ARO-040 検証OK')
