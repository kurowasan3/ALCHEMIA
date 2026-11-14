# POL-RUB-006 検証：陽イオン交換樹脂
import sys
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
assert CalcMolFormula(Chem.MolFromSmiles('C=Cc1ccccc1')) == 'C8H8' and CalcMolFormula(Chem.MolFromSmiles('C=Cc1ccc(C=C)cc1')) == 'C10H10'
unit = 104 * 4 + 130; assert unit == 546
rings = R(546, 10) / unit * 5
so3h = (R(866, 10) - R(546, 10)) / 80
pct = so3h / rings * 100; print(rings, so3h, pct); assert pct == 80
H = R(1, 10) * R(32, 1000); Mg = H / 2; c = Mg / R(20, 1000); print(H, Mg, c); assert c == R(8, 100)
assert H / 2 / R(1, 10) * 1000 == 16
cap = so3h / R(866, 10) * 10; V = cap / 2 / c; print(float(cap), float(V)); assert round(float(V) * 1000) == 289
# R を H に置き換えて式のつり合いを確認
assert check_equation('2HSO3H + Mg^2+ -> (HSO3)2Mg + 2H+') is None
assert not check_item('POL-RUB-006')
print('POL-RUB-006 OK')
