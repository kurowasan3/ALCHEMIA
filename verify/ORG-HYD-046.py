# ORG-HYD-046 検証：C5H10 アルケンの決定と臭素付加物の不斉炭素
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import AllChem
alk = {'1-pentene': 'C=CCCC', '2-pentene': 'CC=CCC', '2-methyl-1-butene': 'C=C(C)CC',
       '3-methyl-1-butene': 'C=CC(C)C', '2-methyl-2-butene': 'CC=C(C)C'}
def h_on_db(s):
    m = Chem.MolFromSmiles(s)
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE:
            return sorted([b.GetBeginAtom().GetTotalNumHs(), b.GetEndAtom().GetTotalNumHs()])
def hydrog(s):
    return Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('=', '')))
def br2(s):
    rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1](Br)[C:2]Br')
    p = rxn.RunReactants((Chem.MolFromSmiles(s),))[0][0]; Chem.SanitizeMol(p)
    return p
for n, s in alk.items():
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m) == 'C5H10'
    p = br2(s)
    cc = Chem.FindMolChiralCenters(p, includeUnassigned=True, useLegacyImplementation=False)
    print(n, 'H on C=C', h_on_db(s), 'H2 ->', hydrog(s), 'Br2 adduct', Chem.MolToSmiles(p), 'stereocenters', len(cc))
# 条件：P: C=C 上の H 合計3・水素付加で枝分かれ → 3-methyl-1-butene
P = [n for n, s in alk.items() if sum(h_on_db(s)) == 3 and hydrog(s) != 'CCCCC']
Q = [n for n, s in alk.items() if h_on_db(s) == [0, 1]]
R = [n for n, s in alk.items() if h_on_db(s) == [0, 2]]
print(P, Q, R)
assert P == ['3-methyl-1-butene'] and Q == ['2-methyl-2-butene'] and R == ['2-methyl-1-butene']
# S: シス-トランス異性体をもつもの
S = []
for n, s in alk.items():
    m = Chem.MolFromSmiles(s)
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
    if len(list(EnumerateStereoisomers(m))) > 1: S.append(n)
print('cis-trans', S); assert S == ['2-pentene']
print('OK')
