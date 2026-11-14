# POL-PRT-018 検証：Ala/Lys ペプチドの組成
import sys
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R, symbols, solve
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import MolWt
from g12lib import check_item

nQ1 = R(487, 1000) / 487; nQ2 = R(974, 1000) / 487
nCO2 = R(224, 1000) / R(224, 10); nN2 = R(784, 10000) / R(224, 10)
n_aa = nCO2 / nQ2; n_N = 2 * nN2 / nQ1
print('アミノ酸数', n_aa, 'N原子数', n_N)
assert n_aa == 5 and n_N == 7
x, y = symbols('x y')
s = solve([x + y - 5, x + 2 * y - 7], [x, y]); print(s)
assert s == {x: 3, y: 2}
assert 89 * 3 + 146 * 2 - 18 * 4 == 487
m = Chem.MolFromSmiles('NC(C)C(=O)NC(C)C(=O)NC(C)C(=O)NC(CCCCN)C(=O)NC(CCCCN)C(=O)O')
print(CalcMolFormula(m), round(MolWt(m), 1), sum(a.GetSymbol() == 'N' for a in m.GetAtoms()))
assert sum(a.GetSymbol() == 'N' for a in m.GetAtoms()) == 7
assert not check_item('POL-PRT-018')
print('POL-PRT-018 OK')
