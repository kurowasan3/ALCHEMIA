# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem2/outh/items/ORG-ARO-042/item.json'
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
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
C = Chem.CanonSmiles


def ozon(smi):
    rw = Chem.RWMol(Chem.MolFromSmiles(smi))
    dbl = [b for b in rw.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE and not b.GetIsAromatic()
           and b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C']
    for b in dbl:
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        rw.RemoveBond(i, j)
        for k in (i, j):
            o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(k, o, Chem.BondType.DOUBLE)
    mm = rw.GetMol(); Chem.SanitizeMol(mm)
    return sorted(C(x) for x in Chem.MolToSmiles(mm).split('.'))


E, F, G = C('O=Cc1ccccc1'), C('O=CCc1ccccc1'), C('CC(=O)C=O')
assert [CalcMolFormula(Chem.MolFromSmiles(x)) for x in (E, F, G)] == ['C7H6O', 'C8H8O', 'C3H4O2']
I = 'c1ccccc1C=CC(C)=CCc1ccccc1'
II = 'c1ccccc1C=C(C)C=CCc1ccccc1'
for s in (I, II):
    assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C18H18'
    assert ozon(s) == sorted([E, F, G])
# G を中央に E・F を両端につなぐ組み立ては 2 通りだけ（G の 2 つのカルボニル炭素のどちらに E がつくか）
assert C(I) != C(II)
nch = {}
for k, s in (('I', I), ('II', II)):
    B = Chem.MolFromSmiles(s.replace('=', ''))
    assert CalcMolFormula(B) == 'C18H22'
    nch[k] = len(Chem.FindMolChiralCenters(B, includeUnassigned=True, useLegacyImplementation=False))
print('B の不斉中心', nch); assert nch == {'I': 0, 'II': 1}
st = {Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(I))}
print('立体異性体', len(st)); assert len(st) == 4
for x in st:
    assert not Chem.FindMolChiralCenters(Chem.MolFromSmiles(x), includeUnassigned=True, useLegacyImplementation=False)
# F の候補（C8H8O の芳香族アルデヒド）と KMnO4 酸化生成物（側鎖→COOH）
cands = {'CC=O側鎖': 'O=CCc1ccccc1', 'o': 'Cc1ccccc1C=O', 'm': 'Cc1cccc(C=O)c1', 'p': 'Cc1ccc(C=O)cc1'}
ox = {'CC=O側鎖': 'OC(=O)c1ccccc1', 'o': 'OC(=O)c1ccccc1C(=O)O', 'm': 'OC(=O)c1cccc(C(=O)O)c1', 'p': 'OC(=O)c1ccc(C(=O)O)cc1'}
benz = C('OC(=O)c1ccccc1')
for k, s in cands.items():
    assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C8H8O'
assert [k for k in ox if C(ox[k]) == benz] == ['CC=O側鎖']
# G：C3H4O2 で CHO と CH3CO- をもつもの
patt_cho = Chem.MolFromSmarts('[CX3H1](=O)[#6]'); patt_ac = Chem.MolFromSmarts('[CH3]C(=O)')
g_c = ['CC(=O)C=O', 'O=CCC=O', 'C=CC(=O)O', 'O=C1CCO1', 'OC=CC=O', 'C=COC=O']
ok = [g for g in g_c if CalcMolFormula(Chem.MolFromSmiles(g)) == 'C3H4O2' and Chem.MolFromSmiles(g).HasSubstructMatch(patt_cho) and Chem.MolFromSmiles(g).HasSubstructMatch(patt_ac)]
assert ok == ['CC(=O)C=O']
print('ORG-ARO-042 検証OK')
