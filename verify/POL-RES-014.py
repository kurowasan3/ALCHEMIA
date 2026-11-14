# POL-RES-014 検証
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
import sympy as sp
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
assert check_equation('C8H12N4 -> 2C4H6N + N2') is None
assert CalcMolFormula(Chem.MolFromSmiles('CC(C)(C#N)N=NC(C)(C)C#N')) == 'C8H12N4'
assert CalcMolFormula(Chem.MolFromSmiles('C[C](C)C#N')) == 'C4H6N'
nN2 = 4.48e-3/22.4; nR = 3.36e-3/14
f = nR/(2*nN2); print('f', f); assert f'{f:.2g}' == '0.6'
M = 1.04*8.31e3*300/(419*0.100); print('M', M); assert f'{M:.3g}' == '6.19e+04'
npol = 10.4/M; print('npol', npol); assert f'{npol:.3g}' == '0.000168'
ratio = nR/npol; print('ratio', ratio)
y = sp.symbols('y'); ysol = sp.solve(sp.Eq(2/(2-y), sp.Rational(240, 168)), y)[0]
print('y', ysol, float(ysol)); assert abs(float(ysol) - 0.60) < 1e-9
ynum = 2 - 2/ratio; print('y(unrounded)', ynum); assert f'{ynum:.2g}' == '0.6'
# 模擬：T 回の停止，再結合 y の割合
T = nR/2; mol = 0.6*T + 2*0.4*T; assert abs(mol - 1.68e-4) < 1e-9
# 連鎖移動を加えると比が下がり y が小さく出る
extra = 0.2e-4
y_app = 2 - 2/(nR/(npol+extra)); print('with transfer y_app', y_app); assert y_app < 0.6
print('OK')
