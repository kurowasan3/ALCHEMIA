import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
from sympy import symbols, solve, log
x = symbols('x', positive=True)
s = solve(x*(x+1e-8) - 1e-14, x)[0]; print(s); approx(float(s), 9.5e-8, 0.005)
h = float(s) + 1e-8; import math; pH = -math.log10(h); print(h, pH); approx(round(pH, 2), 6.98, 1e-9)
approx(round(14-pH, 2), 7.02, 1e-9)
approx(1e-14/1e-5, 1e-9)
assert check_item('THE-ION-012')
