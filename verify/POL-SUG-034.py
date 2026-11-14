# POL-SUG-034 検証：スクロースの部分加水分解物の炭素含有率から加水分解率を求める
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
import sympy as sp
H, C, O = sp.Rational(1), sp.Rational(12), sp.Rational(16)
Msuc = 12 * C + 22 * H + 11 * O; Mhex = 6 * C + 12 * H + 6 * O
cs, ch = 144 / Msuc, 72 / Mhex
print('C% スクロース', float(cs * 100), ' 単糖', float(ch * 100))
assert round(float(cs * 100), 1) == 42.1 and round(float(ch * 100), 1) == 40.0
x = sp.symbols('x')
# 初めのスクロース 1 mol のうち x mol が加水分解 → 質量 342(1-x)+360x，炭素 144 g（不変）
sol = sp.solve(sp.Eq(144 / (Msuc * (1 - x) + 2 * Mhex * x), sp.Rational(410, 1000)), x)[0]
print('加水分解率 x =', float(sol)); assert round(float(sol) * 100) == 51
wsuc = (1 - sol) * Msuc / (Msuc * (1 - sol) + 2 * Mhex * sol)
print('混合物中のスクロースの質量割合', float(wsuc)); assert round(float(wsuc) * 100) == 48
# 同じ加水分解率のマルトースの場合：炭素含有率は同じ式（C12H22O11 → 2C6H12O6）
# Cu2O（還元糖 1 mol → 1 mol とする）：スクロース 2x, マルトース (1-x)+2x = 1+x  (初めの二糖1molあたり)
xs = float(sol)
print('Cu2O mol（二糖1 molあたり）スクロース', 2 * xs, ' マルトース', 1 + xs); assert 1 + xs > 2 * xs
for eq in ['C12H22O11 + H2O -> 2C6H12O6']:
    assert check_equation(eq) is None
print('ALL OK')
