# POL-RES-010 検証：ポリエステル（アジピン酸＋過剰のエチレングリコール）
import sys
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from g12lib import check_item
assert CalcMolFormula(Chem.MolFromSmiles('OC(=O)CCCCC(=O)O')) == 'C6H10O4' and 6*12+10+64 == 146
assert CalcMolFormula(Chem.MolFromSmiles('OCCO')) == 'C2H6O2' and 24+6+32 == 62
na = R(146) / 146; ng = R(682, 10) / 62; print(na, ng); assert ng == R(11, 10)
water = 2 * na * 18; assert water == 36
mols = na + ng - 2 * na; print('分子数', mols); assert mols == R(1, 10)
Mavg = (146 + R(682, 10) - 36) / mols; print('平均分子量', float(Mavg)); assert round(float(Mavg), -2) == 1800
r = na / ng; Xn = (1 + r) / (1 - r); print('Xn', Xn); assert Xn == 21
assert abs(float((na + ng) / mols) - Xn) < 1e-12   # 単量体単位数/分子数 = 数平均重合度
n = R(1) / (2 * 10**4)
dt = R(5) * n / R(1, 10); Pi = n / R(1, 10) * 8300 * 300
print(float(dt), float(Pi)); assert f'{float(dt):.1e}' == '2.5e-03' and f'{float(Pi):.1e}' == '1.2e+03'
assert not check_item('POL-RES-010')
print('POL-RES-010 OK')
