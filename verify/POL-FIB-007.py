# POL-FIB-007 検証：ナイロン610 の分子量
import sys
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R, nsimplify
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item
M = R(282, 1000) * R(83, 10) * 10**3 * 300 / (R(249) * R(1, 10))
print('M =', M, float(M)); assert M == 28200
n = M / 282; print('n =', n); assert n == 100 and 2 * n - 1 == 199
u = Chem.MolFromSmiles('*NCCCCCCNC(=O)CCCCCCCCC(*)=O')
f = CalcMolFormula(u).replace('*2', ''); print(f)
assert f == 'C16H30N2O2' and 16 * 12 + 30 + 2 * 14 + 2 * 16 == 282
assert not check_item('POL-FIB-007')
print('POL-FIB-007 OK')
