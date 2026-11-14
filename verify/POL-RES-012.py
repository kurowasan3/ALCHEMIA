# POL-RES-012 検証：開始剤の効率
import sys
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item
assert CalcMolFormula(Chem.MolFromSmiles('C=C(C)C(=O)OC')) == 'C5H8O2' and 60 + 8 + 32 == 100
assert R(6) / 10 * 100 == 60
assert R(50000) / 100 == 500
chains = R(6) / 50000; eff = chains / R(2, 10000); print(chains, eff); assert eff == R(3, 5)
M = R(6) / (R(4, 10000) * eff); print(M); assert M == 25000
assert not check_item('POL-RES-012')
print('POL-RES-012 OK')
