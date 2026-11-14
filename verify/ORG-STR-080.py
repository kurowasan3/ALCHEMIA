# ORG-STR-080 検証：2-クロロ-4-ヒドロキシ-2-ブテン酸メチルの幾何異性体とラクトン化，E/Z 表示
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import Descriptors
def ez(smi):
    m = Chem.MolFromSmiles(smi)
    from rdkit.Chem import rdCIPLabeler
    rdCIPLabeler.AssignCIPLabels(m)
    labs = [b.GetProp('_CIPCode') for b in m.GetBonds() if b.HasProp('_CIPCode')]
    return CalcMolFormula(m), labs
# 「CH2OH と COOCH3 が同じ側」：C=C の両側で CH2 と C(=O) が cis
# SMILES: OC/C=C(\C(=O)OC)Cl  -> CH2 と C(=O) の関係を RDKit の幾何で確認する
cands = {'X(cis: CH2OH//COOMe)': 'OC/C=C(/Cl)C(=O)OC', 'Y': 'OC/C=C(\\Cl)C(=O)OC'}
from rdkit.Chem import AllChem
import numpy as np
for name, s in cands.items():
    m = Chem.AddHs(Chem.MolFromSmiles(s)); AllChem.EmbedMolecule(m, randomSeed=1)
    conf = m.GetConformer()
    # atoms: 0 O,1 C(H2),2 C3,3 C2,4 Cl,5 C(=O)
    # dihedral C1-C2=C3-C5 (CH2 - C=C - C(=O))
    from rdkit.Chem import rdMolTransforms as T
    dih = T.GetDihedralDeg(conf, 1, 2, 3, 5)
    print(name, s, ez(s), 'dihedral CH2-C=C-C(=O) =', round(dih))
# メチル類似体：Cl の代わりに CH3
for s in ['OC/C=C(/C)C(=O)OC', 'OC/C=C(\\C)C(=O)OC']:
    m = Chem.AddHs(Chem.MolFromSmiles(s)); AllChem.EmbedMolecule(m, randomSeed=1)
    from rdkit.Chem import rdMolTransforms as T
    print('Me analog', s, ez(s), round(T.GetDihedralDeg(m.GetConformer(), 1, 2, 3, 5)))
# ラクトン X'
lac = 'O=C1OCC=C1Cl'
print('lactone', CalcMolFormula(Chem.MolFromSmiles(lac)))
am = {'H': 1.0, 'C': 12.0, 'O': 16.0, 'Cl': 35.5}
def mw(s):
    m = Chem.AddHs(Chem.MolFromSmiles(s)); return sum(am[a.GetSymbol()] for a in m.GetAtoms())
Mx, Ml = mw('OCC=C(Cl)C(=O)OC'), mw(lac)
print('MW ester', Mx, 'MW lactone', Ml)
n = 3.01/Mx; nl = 1.42/Ml
print('n total', n, 'n lactone', nl, 'fraction X', nl/n)
assert abs(nl/n - 0.60) < 0.005
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
print(check_equation('C5H7ClO3 -> C4H3ClO2 + CH3OH'))
