# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem2/outh/items/ORG-HYD-024/item.json'
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
import sympy as sp
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


def products(smi):
    rw = Chem.RWMol(Chem.MolFromSmiles(smi))
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
    for a in rw.GetAtoms():
        a.SetIsAromatic(False); a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    m = rw.GetMol(); Chem.SanitizeMol(m)
    seen, final, stack = set(), set(), [Chem.MolToSmiles(m)]
    while stack:
        s = stack.pop()
        if s in seen:
            continue
        seen.add(s)
        mm = Chem.MolFromSmiles(s)
        three = [r for r in mm.GetRingInfo().AtomRings() if len(r) == 3]
        if not three:
            final.add(s); continue
        bonds = {mm.GetBondBetweenAtoms(r[i], r[(i + 1) % 3]).GetIdx() for r in three for i in range(3)}
        for bi in bonds:
            b = mm.GetBondWithIdx(bi)
            rw = Chem.RWMol(mm); rw.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
            n = rw.GetMol(); Chem.SanitizeMol(n)
            stack.append(Chem.MolToSmiles(n))
    return final


C = Chem.CanonSmiles
cand = {'P': 'C1=CC2C=CC12', 'Q': 'C12C3C1C4C2C34', 'R': 'C1=CC2C3C1C23', 'S': 'C1=CC1C1C=C1'}
want = {
    'P': {C('C1CC2CCC12')},
    'Q': {C('C1CC2CCC12'), C('C1CC2CC1C2')},
    'R': {C('C1CCCCC1'), C('CC1CCCC1'), C('C1CC2CC1C2')},
    'S': {C('CCCCCC'), C('CCCC(C)C'), C('CC(C)C(C)C')},
}
form = {}
for k, s in cand.items():
    assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C6H6'
    P = {C(p) for p in products(s)}
    form[k] = sorted({CalcMolFormula(Chem.MolFromSmiles(p)) for p in P})
    print(k, sorted(P), form[k])
    assert P == want[k], k
assert form == {'P': ['C6H10'], 'Q': ['C6H10'], 'R': ['C6H10', 'C6H12'], 'S': ['C6H14']}
# X：C6H12 2 種類（シクロヘキサンを含む）＋C6H10 1 種類 → R のみ
for k in cand:
    P = {C(p) for p in products(cand[k])}
    f = [CalcMolFormula(Chem.MolFromSmiles(p)) for p in P]
    isX = f.count('C6H12') == 2 and f.count('C6H10') == 1 and C('C1CCCCC1') in P
    isY = len(P) == 3 and all(x == 'C6H14' for x in f)
    print(k, 'X候補' if isX else '', 'Y候補' if isY else '')
    assert isX == (k == 'R') and isY == (k == 'S')
# 2 つの三員環が結合を共有するのは R だけ
for k, s in cand.items():
    rings = [set(r) for r in Chem.GetSymmSSSR(Chem.MolFromSmiles(s)) if len(r) == 3]
    share = any(len(a & b) == 2 for i, a in enumerate(rings) for b in rings[i + 1:])
    assert share == (k == 'R')
n = sp.symbols('n')
sol = sp.solve(sp.Eq(2 * n + 3 * (sp.Rational(1, 100) - n), sp.Rational(28, 1000)), n)
print('C6H10', sol); assert sol == [sp.Rational(2, 1000)]
print('ORG-HYD-024 検証OK')
