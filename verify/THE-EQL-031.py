import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
from sympy import symbols, solve, sqrt
Kp = 1.6e9
p = (Kp)**0.5; approx(2*p, 8.0e4)
x = symbols('x'); r = sorted(solve(x*(1.0e5-x) - Kp, x)); print(r); approx(float(r[0]), 2.0e4); approx(float(r[1]), 8.0e4)
assert check_item('THE-EQL-031')
