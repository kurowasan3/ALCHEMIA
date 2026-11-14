# POL-RES-011 検証：ポリエチレンの結晶化度
import sys
sys.path.insert(0, '/home/claude/chem/work/g12')
from sympy import Rational as R, symbols, solve
from g12lib import check_item
x = symbols('x')
res = {}
for name, dd in (('P', R(92, 100)), ('Q', R(95, 100))):
    s = solve(1 / dd - (x / 1 + (1 - x) / R(85, 100)), x)[0]
    res[name] = float(s); print(name, float(s))
assert round(res['P'] * 100) == 51 and round(res['Q'] * 100) == 70
assert not check_item('POL-RES-011')
print('POL-RES-011 OK')
