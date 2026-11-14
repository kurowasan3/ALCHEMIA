"""THE-MOL-045 検証"""
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
n = 1000 / 2.0
VA = n * 8.31e3 * 298 / 3.5e7
VB = 1.0 / 0.071
mC = n / 5 * (10 * 12 + 18 * 1.0) / 1000
VC = mC / 0.88
mD = n * (24.3 + 2.0) / 1000
VD = mD / 1.45
print(VA, VB, mC, VC, mD, VD)
r = lambda x: float(f'{x:.2g}')
assert (r(VA), r(VB), r(mC), r(VC), r(VD)) == (35, 14, 14, 16, 9.1)
assert sorted([('A', VA), ('B', VB), ('C', VC), ('D', VD)], key=lambda t: t[1])[0][0] == 'D'
assert [k for k, v in sorted([('A', VA), ('B', VB), ('C', VC), ('D', VD)], key=lambda t: t[1])] == ['D', 'B', 'C', 'A']
assert max([('A', 1.0), ('B', 1.0), ('C', mC), ('D', mD)], key=lambda t: t[1])[0] == 'C'
assert check_equation('C10H8 + 5H2 -> C10H18') is None
assert CalcMolFormula(Chem.MolFromSmiles('c1ccc2ccccc2c1')) == 'C10H8'
assert CalcMolFormula(Chem.MolFromSmiles('C1CCC2CCCCC2C1')) == 'C10H18'
print('THE-MOL-045 OK', 2.0 / (mC) * 100 / 1000 * 500, 1.0 / mD * 100)
