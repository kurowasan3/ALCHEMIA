# ORG-HYD-014 検証：アルキンへの水の付加（マルコフニコフ則＋エノールの異性化）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def canon(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s))
def hydrate_alkyne(smi):
    """三重結合 1 本に H2O を付加：OH は H の少ない炭素へ（等しければ両方），エノール→ケト"""
    m = Chem.MolFromSmiles(smi)
    b = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.TRIPLE][0]
    a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
    h1, h2 = a1.GetTotalNumHs(), a2.GetTotalNumHs()
    targets = [a1] if h1 < h2 else [a2] if h2 < h1 else [a1, a2]
    prods = set()
    for t in targets:
        rw = Chem.RWMol(m)
        rw.GetBondBetweenAtoms(a1.GetIdx(), a2.GetIdx()).SetBondType(Chem.BondType.SINGLE)
        o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(t.GetIdx(), o, Chem.BondType.DOUBLE)   # ケト形（エノールの異性化後）
        for a in rw.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
        mm = rw.GetMol(); Chem.SanitizeMol(mm); prods.add(Chem.MolToSmiles(mm))
    return sorted(prods)
res = {k: hydrate_alkyne(s) for k, s in {'acetylene': 'C#C', 'propyne': 'C#CC', '1-butyne': 'C#CCC', '2-butyne': 'CC#CC'}.items()}
print(res)
assert res['acetylene'] == [canon('CC=O')] and res['propyne'] == [canon('CC(C)=O')]
assert res['1-butyne'] == res['2-butyne'] == [canon('CCC(C)=O')]
ald = Chem.MolFromSmarts('[CX3H1](=O)')
assert [k for k, v in res.items() if any(Chem.MolFromSmiles(x).HasSubstructMatch(ald) for x in v)] == ['acetylene']
# エノール中間体の分子式＝ケト形と同じ（異性体）
assert CalcMolFormula(Chem.MolFromSmiles('C=C(O)C')) == CalcMolFormula(Chem.MolFromSmiles('CC(C)=O')) == 'C3H6O'
for eq in ['CH#CH + H2O -> CH3CHO', 'CH#CCH3 + H2O -> CH3COCH3', 'CH3COOCH=CH2 + H2O -> CH3COOH + CH3CHO',
           'HOCH2CH2OH -> CH3CHO + H2O', 'CH3CH2C#CH + H2O -> CH3CH2COCH3']:
    assert check_equation(eq) is None, eq
print('ORG-HYD-014 OK')
