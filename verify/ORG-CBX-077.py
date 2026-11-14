# ORG-CBX-077 検証：メタリルアルコール -> (Br2) -> (KMnO4) -> (Zn) メタクリル酸
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import FindMolChiralCenters
mols = {'A': 'C=C(C)CO', 'B': 'BrCC(C)(Br)CO', 'C': 'BrCC(C)(Br)C(=O)O', 'P': 'C=C(C)C(=O)O'}
am = {'H': 1.0, 'C': 12.0, 'O': 16.0, 'Br': 79.9}
for k, s in mols.items():
    m = Chem.MolFromSmiles(s)
    mh = Chem.AddHs(m)
    mw = sum(am[a.GetSymbol()] for a in mh.GetAtoms())
    cc = FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
    print(k, CalcMolFormula(m), round(mw, 1), 'chiral centers', cc)
for eq in ['CH2=C(CH3)CH2OH + Br2 -> CH2BrCBr(CH3)CH2OH',
           'CH2BrCBr(CH3)COOH + Zn -> CH2=C(CH3)COOH + ZnBr2',
           '5CH2BrCBr(CH3)CH2OH + 4KMnO4 + 6H2SO4 -> 5CH2BrCBr(CH3)COOH + 4MnSO4 + 2K2SO4 + 11H2O']:
    print(eq, '->', check_equation(eq))
    assert check_equation(eq) is None
n = 3.60/72.0
y = 0.85*0.70*0.90
m = n*y*86.0
print('overall yield', y, 'mass', m)
assert f'{m:.2f}' == '2.30'
