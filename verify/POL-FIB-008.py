# POL-FIB-008 検証：ナイロン66 の分子量（浸透圧と末端基）
import sys, re
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R, symbols, solve
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item
f = re.sub(r'\*\d*', '', CalcMolFormula(Chem.MolFromSmiles('*NCCCCCCNC(=O)CCCCC(*)=O'))); print(f)
assert f == 'C12H22N2O2' and 12*12 + 22 + 28 + 32 == 226
M1 = R(113, 100) * 8300 * 300 / (R(2490) * R(1, 10)); print(M1); assert M1 == 11300
nNH2 = R(1, 100) * R(7, 1000); M2 = R(565, 1000) / nNH2; print(float(M2)); assert round(float(M2), -2) == 8100
nmol = R(565, 1000) / M1; x = symbols('x'); sx = solve(nmol * (2 * x + 1 - x) - nNH2, x)[0]; print(nmol, sx); assert sx == R(2, 5)
assert M1 / 226 == 50
assert not check_item('POL-FIB-008')
print('POL-FIB-008 OK')
