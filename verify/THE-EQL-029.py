import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
from sympy import symbols, solve
approx(0.010/20*2/0.50*1000, 2.0)       # 空試験 2.0 mL
ac = 0.50*(27.0-2.0)/1000*20; approx(ac, 0.25)
x = 1-ac; K = x*x/(ac*(2.00-x)); approx(K, 1.8, 1e-6)
ac2 = 0.50*27.0/1000*20; x2 = 1-ac2; K2 = x2*x2/(ac2*(2.00-x2)); print(K2); approx(K2, 1.55, 0.005); assert K2 < K
a = symbols('a'); s = solve(0.81/(0.10*(a-0.90)) - 1.8, a)[0]; approx(float(s), 5.4)
assert check_item('THE-EQL-029')
