# POL-SUG-052 検証：マルトース・スクロース・ラクトース混合物の加水分解生成物の比と Cu2O
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, Rational as R, N
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
a, b, c = symbols('a b c', positive=True)
sol = solve([ (2*a + b + c) - 9*c/3, b - 2*c/3 ], [a, b], dict=True)[0]
print(sol)   # c=3 とおくと a=2, b=2
ratio = [sol[a].subs(c, 3), sol[b].subs(c, 3), 3]
assert ratio == [2, 2, 3]
Mdi = 12*12 + 22*R(1) + 11*16; Mcu2o = 2*64 + 16
assert Mdi == 342 and Mcu2o == 144
n = R('4.79')/Mdi; print('n total =', N(n))
na, nb, nc = n*R(2, 7), n*R(2, 7), n*R(3, 7)
before = (na + nc)*Mcu2o; after = 2*n*Mcu2o
print('Cu2O before =', N(before), 'after =', N(after))
assert round(float(before), 2) == 1.44 and round(float(after), 2) == 4.03
# (4) アルドース:ケトース比と加水分解前の Cu2O はどちらも a+c だけに依存
for (x, z) in [(2, 3), (4, 1), (0, 5)]:
    y = 2
    aldo, keto, red = 2*x + y + 2*z, y, x + z
    print((x, y, z), 'aldose:ketose =', aldo, ':', keto, 'reducing =', red)
# 分子式の確認
smi = {'maltose': 'OC[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O',
       'sucrose': 'OC[C@H]1O[C@H](O[C@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
       'lactose': 'OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@H]1O'}
for k, s in smi.items():
    assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C12H22O11', k
for e in ['C12H22O11 + H2O -> 2C6H12O6']:
    assert check_equation(e) is None
print('OK')
