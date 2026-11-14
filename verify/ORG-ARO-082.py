# ORG-ARO-082 検証：環状イミドの開環生成物の数と構造決定
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
def openings(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    for b in m.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if not b.IsInRing(): continue
        pair = {a1.GetSymbol(), a2.GetSymbol()}
        if pair != {'C', 'N'}: continue
        c = a1 if a1.GetSymbol() == 'C' else a2
        if not any(n.GetSymbol() == 'O' and m.GetBondBetweenAtoms(c.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2 for n in c.GetNeighbors()): continue
        rw = Chem.RWMol(m); rw.RemoveBond(a1.GetIdx(), a2.GetIdx())
        o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(c.GetIdx(), o, Chem.BondType.SINGLE)
        p = rw.GetMol(); Chem.SanitizeMol(p); out.add(Chem.MolToSmiles(p))
    return out
fig = {'a': ('CN1C(=O)c2ccccc2C1=O', 1), 'b': ('CN1C(=O)c2cccc(C)c2C1=O', 2), 'c': ('CN1C(=O)c2ccc(C)cc2C1=O', 2),
       'd': ('CN1C(=O)CCC1=O', 1), 'e': ('CN1C(=O)CC(C)C1=O', 2)}
for k, (s, n) in fig.items():
    o = openings(s); print(k, F(Chem.MolFromSmiles(s)), len(o)); assert len(o) == n
same = [k for k, (s, n) in fig.items() if F(Chem.MolFromSmiles(s)) == 'C10H9NO2']; assert same == ['b', 'c']
Z = 'CN1C(=O)Cc2ccccc2C1=O'; assert F(Chem.MolFromSmiles(Z)) == 'C10H9NO2'
W = openings(Z); print(W); assert len(W) == 2 and all(F(Chem.MolFromSmiles(w)) == 'C10H11NO3' for w in W)
ring = [len(r) for r in Chem.GetSymmSSSR(Chem.MolFromSmiles(Z))]; assert sorted(ring) == [6, 6]
assert F(Chem.MolFromSmiles('OC(=O)Cc1ccccc1C(=O)O')) == 'C9H8O4'
for e in ['C10H9NO2 + H2O -> C10H11NO3', 'C10H11NO3 + H2O -> C9H8O4 + CH5N']:
    assert check_equation(e) is None
# C9H8O4 の候補と酸化生成物（ベンジル位に H をもつ側鎖 → COOH）
cands = {'3-メチルフタル酸': ('Cc1cccc(C(=O)O)c1C(=O)O', 'C9H6O6'), '4-メチルフタル酸': ('Cc1ccc(C(=O)O)c(C(=O)O)c1', 'C9H6O6'),
         'フェニルマロン酸': ('OC(=O)C(C(=O)O)c1ccccc1', 'C7H6O2'), 'ホモフタル酸': ('OC(=O)Cc1ccccc1C(=O)O', 'C8H6O4')}
for k, (s, ox) in cands.items(): assert F(Chem.MolFromSmiles(s)) == 'C9H8O4'
print('OK')
