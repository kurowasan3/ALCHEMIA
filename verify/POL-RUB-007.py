# POL-RUB-007 検証：クロロプレンゴム
import sys, re
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item
F = lambda s: re.sub(r'\*\d*', '', CalcMolFormula(Chem.MolFromSmiles(s)))
assert F('C=CC#C') == 'C4H4' and F('C=CC(Cl)=C') == 'C4H5Cl' and F('*CC(Cl)=CC*') == 'C4H5Cl'
unit = 12 * 4 + R(10, 10) * 5 + R(355, 10); assert unit == R(885, 10)
V = R(177) / unit * 2 * R(224, 10); print(V); assert V == R(896, 10)
# 1,2-付加：主鎖（* と * を結ぶ原子）に二重結合がない／1,4-付加：ある
def main_chain_double(s):
    m = Chem.MolFromSmiles(s)
    stars = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == '*']
    path = Chem.GetShortestPath(m, *stars)
    return any(m.GetBondBetweenAtoms(path[i], path[i + 1]).GetBondTypeAsDouble() == 2 for i in range(len(path) - 1))
print(main_chain_double('*CC(*)C=C'), main_chain_double('*CC=CC*'))
assert not main_chain_double('*CC(*)C=C') and main_chain_double('*CC=CC*')
assert not check_item('POL-RUB-007')
print('POL-RUB-007 OK')
