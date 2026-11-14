# ORG-ALC-033 検証：アセトニド保護による選択的酸化
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
for e in ['C4H10O3 + C3H6O -> C7H14O3 + H2O', 'C7H12O4 + H2O -> C4H8O4 + C3H6O', 'C4H8O4 -> C4H6O3 + H2O']:
    assert check_equation(e) is None
D = {'OCCC(O)CO': ('C4H10O3', []), 'CC1(C)OCC(CCO)O1': ('C7H14O3', [5]), 'CC1(C)OCCC(CO)O1': ('C7H14O3', [6]),
     'CC1(C)OCC(CC(=O)O)O1': ('C7H12O4', [5]), 'OCC(O)CC(=O)O': ('C4H8O4', []), 'O=C1CC(O)CO1': ('C4H6O3', [5]),
     'O=C1CC(CO)O1': ('C4H6O3', [4])}
for s, (f, r) in D.items():
    m = Chem.MolFromSmiles(s); assert F(m) == f and [len(x) for x in Chem.GetSymmSSSR(m)] == r, s
print('OK')
