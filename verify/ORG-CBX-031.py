# 検証スクリプト（グループ g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-CBX-031/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)

# ---- 計算
from sympy import Rational as R, symbols, solve

# ---- 計算
from sympy import Rational as R, nsimplify
mC = R(440,10)*12/44; mH = R(108,10)*2/18; mO = R(260,10) - mC - mH
r = [mC/12, mH/1, mO/16]; g = r[0]/5
print('C:H:O =', [x/g for x in r])
assert [x/g for x in r] == [5, 6, 4]
nNaOH = R(1,10)*R(20,1000)
cands = [n for n in range(1, 10) if (R(130,1000)/(nNaOH/n)) % 130 == 0 and R(130,1000)/(nNaOH/n) <= 200]
print('価数の候補', cands); assert cands == [2]
nX = R(130,1000)/130; nBr2 = R(160,1000)/160; assert nX == nBr2

# ---- RDKit：C5H6O4 の鎖式ジカルボン酸（C=C 1個）を網羅
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as Dm
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from itertools import combinations_with_replacement as cwr
base = Chem.MolFromSmiles('C=CC'); S = set()
for i, j in cwr(range(3), 2):
    rw = Chem.RWMol(base)
    for k in (i, j):
        c = rw.AddAtom(Chem.Atom(6)); o1 = rw.AddAtom(Chem.Atom(8)); o2 = rw.AddAtom(Chem.Atom(8))
        rw.AddBond(k, c, Chem.BondType.SINGLE); rw.AddBond(c, o1, Chem.BondType.DOUBLE); rw.AddBond(c, o2, Chem.BondType.SINGLE)
    try:
        m = rw.GetMol(); Chem.SanitizeMol(m); S.add(Chem.MolToSmiles(m))
    except Exception:
        pass
assert all(Dm.CalcMolFormula(Chem.MolFromSmiles(s)) == 'C5H6O4' for s in S)
print('構造異性体', len(S)); assert len(S) == 5
def hyd(s):
    rw = Chem.RWMol(Chem.MolFromSmiles(s))
    for b in rw.GetBonds():
        if b.GetBondTypeAsDouble() == 2 and b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C':
            b.SetBondType(Chem.BondType.SINGLE)
    m = rw.GetMol(); Chem.SanitizeMol(m); return m
def nch(m): return len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
ok = []
for s in S:
    m = Chem.MolFromSmiles(s)
    ct = len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m)}) == 2
    h = hyd(s)
    print(s, 'シストランス', ct, '付加体', Chem.MolToSmiles(h), '不斉', nch(h))
    if ct and nch(h) == 1: ok.append(s)
assert ok == [Chem.CanonSmiles('CC(=CC(=O)O)C(=O)O')], ok
same = [s for s in S if s != ok[0] and Chem.MolToSmiles(hyd(s)) == Chem.MolToSmiles(hyd(ok[0]))]
assert same == [Chem.CanonSmiles('C=C(CC(=O)O)C(=O)O')], same
# シス形（2つの COOH が同じ側）は Z
from rdkit.Chem import rdCIPLabeler
m = Chem.MolFromSmiles('OC(=O)/C(C)=C\\C(=O)O'); rdCIPLabeler.AssignCIPLabels(m)
assert any(b.GetPropsAsDict().get('_CIPCode') == 'Z' for b in m.GetBonds())
print('OK')
