# ORG-STR-011 検証：ヒドロキシ酸 B と A の構造
import sys, itertools
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item
mC = R(440, 10) * 12 / 44; mH = R(180, 10) * 2 / 18; mO = 26 - mC - mH
r = (mC / 12, mH / 1, mO / 16); k = r[0] / 4
print(mC, mH, mO, [x / k for x in r]); assert [x / k for x in r] == [4, 8, 3]
assert 4 * 12 + 8 + 3 * 16 == 104
# C4 骨格（直鎖・枝分かれ）の各炭素に COOH 炭素と OH を置く総当たり
cands = set()
for skel in ('CCCC', 'CC(C)C'):
    m = Chem.MolFromSmiles(skel)
    n = m.GetNumAtoms()
    for i in range(n):               # i 番目の炭素をカルボキシ基の炭素にする（末端のみ可）
        if m.GetAtomWithIdx(i).GetDegree() != 1: continue
        for j in range(n):
            if j == i: continue
            rw = Chem.RWMol(m)
            o1 = rw.AddAtom(Chem.Atom(8)); o2 = rw.AddAtom(Chem.Atom(8)); o3 = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(i, o1, Chem.BondType.DOUBLE); rw.AddBond(i, o2, Chem.BondType.SINGLE); rw.AddBond(j, o3, Chem.BondType.SINGLE)
            try:
                Chem.SanitizeMol(rw)
            except Exception:
                continue
            cands.add(Chem.MolToSmiles(rw))
cands = {c for c in cands if CalcMolFormula(Chem.MolFromSmiles(c)) == 'C4H8O3'}
iodo = Chem.MolFromSmarts('[CH3][CH1]([OX2H1])')
res = []
for c in sorted(cands):
    m = Chem.MolFromSmiles(c)
    ch = bool(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
    io = m.HasSubstructMatch(iodo)
    res.append((c, ch, io)); print(c, '不斉' if ch else '', 'ヨードホルム+' if io else '')
assert len(cands) == 5
B = [c for c, ch, io in res if ch and io]; assert B == [Chem.CanonSmiles('CC(O)CC(=O)O')]
q4 = sorted(c for c, ch, io in res if ch and not io)
assert q4 == sorted([Chem.CanonSmiles('CCC(O)C(=O)O'), Chem.CanonSmiles('OCC(C)C(=O)O')])
A = Chem.MolFromSmiles('CC(O)CC(=O)NCC(=O)O'); print(CalcMolFormula(A))
assert CalcMolFormula(A) == 'C6H11NO4' and 6 * 12 + 11 + 14 + 4 * 16 == 161
assert not check_item('ORG-STR-011')
print('ORG-STR-011 OK')
