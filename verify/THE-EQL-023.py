import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
from sympy import symbols, solve, Rational as Q
x = symbols('x', positive=True)
K = 36
print('v1/v2 at HI 0.5', K/(Q(1,2)**2/Q(1,4)**2)); approx(float(K/(Q(1,2)**2/Q(1,4)**2)), 9.0)
s = [r for r in solve((2*x)**2 - K*(Q(1,2)-x)**2, x) if r < Q(1,2)][0]; print('x', s); approx(float(2*s), 0.75)
q = Q(8,10)**2/(Q(1,10)**2); approx(float(K/q), 0.5625)
y = [r for r in solve((Q(8,10)-2*x)**2 - K*(Q(1,10)+x)**2, x) if r < Q(4,10)][0]; print('y', y); approx(float(Q(8,10)-2*y), 0.75)
assert check_item('THE-EQL-023')
