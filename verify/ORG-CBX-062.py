# ORG-CBX-062 検証：燃焼分析 → C7H12O4，ジカルボン酸モノエステルの網羅と条件による絞り込み
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import MolWt
from sympy import Rational as R, nsimplify

# 1) 燃焼分析
mA, mCO2, mH2O = R('32.0'), R('61.6'), R('21.6')
mC = mCO2*12/44; mH = mH2O*2/18; mO = mA - mC - mH
nC, nH, nO = mC/12, mH/1, mO/16
print('mC,mH,mO =', float(mC), float(mH), float(mO))
r = [nC/nO, nH/nO, nO/nO]
print('C:H:O =', [float(x*4) for x in r])
assert [x*4 for x in r] == [7, 12, 4]
assert 7*12+12+4*16 == 160   # 組成式量 160，分子量 200 以下なので n=1

# 2) 網羅：飽和鎖式の骨格に COOH を 2 個（ジカルボン酸），OH を 1 個（アルコール）
skel = ['C', 'CC', 'CCC', 'CCCC', 'CC(C)C', 'CCCCC', 'CC(C)CC', 'CC(C)(C)C']
def subst(smi, grp, k):
    m = Chem.MolFromSmiles(smi); out = set()
    idx = [a.GetIdx() for a in m.GetAtoms()]
    import itertools
    for comb in itertools.combinations_with_replacement(idx, k):
        rw = Chem.RWMol(m); ok = True
        for i in comb:
            a = rw.GetAtomWithIdx(i)
            if a.GetTotalNumHs() == 0: ok = False; break
            n = rw.AddAtom(Chem.Atom(6 if grp=='COOH' else 8))
            rw.AddBond(i, n, Chem.BondType.SINGLE)
            if grp == 'COOH':
                o1 = rw.AddAtom(Chem.Atom(8)); o2 = rw.AddAtom(Chem.Atom(8))
                rw.AddBond(n, o1, Chem.BondType.DOUBLE); rw.AddBond(n, o2, Chem.BondType.SINGLE)
            Chem.SanitizeMol(rw)
        if ok:
            try:
                Chem.SanitizeMol(rw); out.add(Chem.MolToSmiles(rw))
            except Exception: pass
    return out
diacids = {'OC(=O)C(=O)O'}
for s in skel: diacids |= subst(s, 'COOH', 2)
alcohols = set()
for s in skel: alcohols |= subst(s, 'OH', 1)
def chiral(smi):
    return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(smi), includeUnassigned=True, useLegacyImplementation=False)) > 0
patt_cooh = Chem.MolFromSmarts('C(=O)[OX2H1]')
cands = set()
rxn_smarts = '[C:1](=[O:2])[OX2H1:3]'
for d in diacids:
    md = Chem.MolFromSmiles(d)
    for match in md.GetSubstructMatches(patt_cooh):
        for al in alcohols:
            ma = Chem.MolFromSmiles(al)
            # エステル化：d の match[2] の O を削り，アルコールの O と結合
            combo = Chem.RWMol(Chem.CombineMols(md, ma))
            off = md.GetNumAtoms()
            oh = [a.GetIdx() for a in ma.GetAtoms() if a.GetSymbol()=='O'][0] + off
            combo.AddBond(match[0], oh, Chem.BondType.SINGLE)
            combo.RemoveAtom(match[2])
            m2 = combo.GetMol(); Chem.SanitizeMol(m2)
            if CalcMolFormula(m2) == 'C7H12O4':
                cands.add((Chem.MolToSmiles(m2), d, al))
print('モノエステル C7H12O4 の候補数（酸・アルコールの組）:', len(cands))
iodo = Chem.MolFromSmarts('[CH3][CX4;H1,H2][OX2H1]')
prim = Chem.MolFromSmarts('[CX4;H2,H3][OX2H1]')
ok = []
for A, B, Cc in sorted(cands):
    mc = Chem.MolFromSmiles(Cc)
    conds = [chiral(A), not chiral(B), not chiral(Cc), mc.HasSubstructMatch(iodo)]
    if all(conds):
        ok.append((A, B, Cc, mc.HasSubstructMatch(prim)))
for x in ok: print('不斉A・アキラルB,C・ヨードホルム陽性:', x)
final = [x for x in ok if x[3]]
print('さらに C が第一級:', final)
assert len(final) == 1
assert Chem.MolToSmiles(Chem.MolFromSmiles(final[0][0])) == Chem.MolToSmiles(Chem.MolFromSmiles('CCOC(=O)C(CC)C(=O)O'))
sec = [x for x in ok if not x[3]]
assert len(sec) == 1 and Chem.MolToSmiles(Chem.MolFromSmiles(sec[0][0])) == Chem.MolToSmiles(Chem.MolFromSmiles('CC(C)OC(=O)C(C)C(=O)O'))
print('MolWt A =', MolWt(Chem.MolFromSmiles('CCOC(=O)C(CC)C(=O)O')))
for e in ['C7H12O4 + 8O2 -> 7CO2 + 6H2O',
          'C7H12O4 + H2O -> C5H8O4 + C2H6O']:
    assert check_equation(e) is None, e
print('OK')
