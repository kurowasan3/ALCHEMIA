ITEM = '/home/claude/chemk/outh/items/POL-SUG-056/item.json'
# 検証スクリプト（グループ g04）  python3 -I この.py
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
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
F = lambda s: CalcMolFormula(Chem.MolFromSmiles(s)); can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
nC = lambda s: len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True))
DHA = 'OCC(=O)CO'
fru = 'OCC(=O)C(O)C(O)C(O)CO'; c4 = 'OCC(=O)C(O)CO'; c5 = 'OCC(=O)C(O)C(O)CO'
assert F(fru) == 'C6H12O6' and nC(fru) == 3 and nC(c4) == 1 and nC(c5) == 2
assert F(c4) == 'C4H8O4' and F(c5) == 'C5H10O5'
# (2) C5 を固定した立体異性体
fixed = Chem.MolFromSmiles('OCC(=O)C(O)C(O)[C@H](O)CO')
iso = set(Chem.MolToSmiles(m) for m in EnumerateStereoisomers(fixed, options=StereoEnumerationOptions(onlyUnassigned=True)))
print(len(iso), iso); assert len(iso) == 4
# 逆アルドールの切断：C=O 炭素の隣の炭素 α と，その先で OH をもつ β の結合を切る
def retro(s):
    m = Chem.MolFromSmiles(s); res = set()
    for a in m.GetAtoms():
        if a.GetSymbol() != 'C': continue
        if not any(nb.GetSymbol() == 'O' and m.GetBondBetweenAtoms(a.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 2 for nb in a.GetNeighbors()): continue
        for al in a.GetNeighbors():
            if al.GetSymbol() != 'C': continue
            for be in al.GetNeighbors():
                if be.GetIdx() == a.GetIdx() or be.GetSymbol() != 'C': continue
                ohs = [x for x in be.GetNeighbors() if x.GetSymbol() == 'O' and x.GetTotalNumHs() == 1]
                if not ohs: continue
                rw = Chem.RWMol(m)
                rw.RemoveBond(al.GetIdx(), be.GetIdx())
                rw.GetBondBetweenAtoms(be.GetIdx(), ohs[0].GetIdx()).SetBondType(Chem.BondType.DOUBLE)
                for at in rw.GetAtoms(): at.SetNoImplicit(False); at.SetNumExplicitHs(0)
                mm = rw.GetMol(); Chem.SanitizeMol(mm)
                res.add(tuple(sorted(Chem.MolToSmiles(f) for f in Chem.GetMolFrags(mm, asMols=True))))
    return res
sedo = 'OCC(=O)C(O)C(O)C(O)C(O)CO'; glc = 'O=CC(O)C(O)C(O)C(O)CO'
r1 = retro(sedo); r2 = retro(glc); r3 = retro(fru); print(r1, r2, r3)
ery = can('O=CC(O)C(O)CO')
assert r1 == {tuple(sorted([can(DHA), ery]))}
assert r2 == {tuple(sorted([can('O=CCO'), ery]))}
assert r3 == {tuple(sorted([can(DHA), can('O=CC(O)CO')]))}
assert F(sedo) == 'C7H14O7'

print('全検証OK:', ITEM)
