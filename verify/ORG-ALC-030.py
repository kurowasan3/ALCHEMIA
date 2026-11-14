# ORG-ALC-030 検証：C4H10O2 の異性体（ジオール・エーテルアルコール）の列挙と条件による構造決定
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import FindMolChiralCenters
from itertools import combinations
skel = {'butane': 'CCCC', 'isobutane': 'CC(C)C'}
cands = set()
# ジオール：異なる2つの炭素に OH
for s in skel.values():
    m = Chem.MolFromSmiles(s)
    for i, j in combinations(range(m.GetNumAtoms()), 2):
        rw = Chem.RWMol(m)
        for k in (i, j):
            o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(k, o, Chem.BondType.SINGLE)
        cands.add(Chem.MolToSmiles(rw))
# エーテルアルコール：R-O-R'（炭素合計4）に OH を1つ（エーテル O と同じ炭素は除く）
alkyl = {1: ['C'], 2: ['CC'], 3: ['CCC', 'C(C)C']}
for n1 in (1, 2):
    for a in alkyl[n1]:
        for b in alkyl[4 - n1]:
            smi = a + 'O' + b
            m = Chem.MolFromSmiles(smi)
            oi = [x.GetIdx() for x in m.GetAtoms() if x.GetSymbol() == 'O'][0]
            nb = [x.GetIdx() for x in m.GetAtomWithIdx(oi).GetNeighbors()]
            for c in range(m.GetNumAtoms()):
                if m.GetAtomWithIdx(c).GetSymbol() != 'C' or c in nb: continue
                rw = Chem.RWMol(m); o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(c, o, Chem.BondType.SINGLE)
                cands.add(Chem.MolToSmiles(rw))
cands = sorted(cands)
for c in cands: assert CalcMolFormula(Chem.MolFromSmiles(c)) == 'C4H10O2'
print(len(cands), cands)
def nOH(s): return len(Chem.MolFromSmiles(s).GetSubstructMatches(Chem.MolFromSmarts('[CX4][OX2H]')))
def iodoform(s): return Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts('[CH3][CX4H1][OX2H]'))
def oxidize(s):
    m = Chem.RWMol(Chem.MolFromSmiles(s)); Chem.Kekulize(m)
    prim = m.GetSubstructMatches(Chem.MolFromSmarts('[CX4;H2,H3][OX2H]'))      # 第一級 → COOH
    sec = [(t[0], t[3]) for t in m.GetSubstructMatches(Chem.MolFromSmarts('[CX4;H1]([#6])([#6])[OX2H]'))]  # 第二級 → C=O
    for c, o in prim:
        oo = m.AddAtom(Chem.Atom(8)); m.AddBond(c, oo, Chem.BondType.DOUBLE)
    for c, o in sec:
        m.GetBondBetweenAtoms(c, o).SetBondType(Chem.BondType.DOUBLE)
    mm = m.GetMol(); 
    for a in mm.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    Chem.SanitizeMol(mm); return Chem.MolToSmiles(mm)
def branched(s):
    m = Chem.MolFromSmiles(s)
    return any(sum(1 for n in a.GetNeighbors() if n.GetSymbol() == 'C') > 2 for a in m.GetAtoms() if a.GetSymbol() == 'C')
def chiral(s): return len(FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False)) > 0
A = [s for s in cands if nOH(s) == 2 and not iodoform(s) and not branched(s) and Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts('[CH2]O')) and CalcMolFormula(Chem.MolFromSmiles(oxidize(s))) == 'C4H6O4']
B = [s for s in cands if nOH(s) == 2 and iodoform(s) and CalcMolFormula(Chem.MolFromSmiles(oxidize(s))) == 'C4H6O3']
C = [s for s in cands if nOH(s) == 1 and chiral(s) and CalcMolFormula(Chem.MolFromSmiles(oxidize(s))) == 'C4H8O3' and chiral(oxidize(s))]
print('A', A, [oxidize(s) for s in A]); print('B', B, [oxidize(s) for s in B]); print('C', C, [oxidize(s) for s in C])
assert A == ['OCCCCO'] and B == ['CC(O)CCO'] and len(C) == 1 and Chem.CanonSmiles(C[0]) == Chem.CanonSmiles('COC(C)CO')
# 候補の分類（参考）：ジオールで酸化生成物 C4H6O4 になるもの（枝分かれ不問）
print('C4H6O4 diols', [s for s in cands if nOH(s) == 2 and CalcMolFormula(Chem.MolFromSmiles(oxidize(s))) == 'C4H6O4'])
print('B has chiral', chiral('CC(O)CCO'))
