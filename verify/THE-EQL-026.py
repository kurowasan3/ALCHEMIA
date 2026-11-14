import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
from sympy import symbols, simplify
approx((0.020+0.020)*5.0, 0.20)
Kc = 0.020**2*0.010/0.020**2; approx(Kc, 0.010); approx(Kc*8.3e3*1000, 8.3e4)
R, T, a, b, c = symbols('R T a b c', positive=True)
p = lambda s: s*R*T
print(simplify((p(c)**2/(p(a)*p(b))) / (c**2/(a*b))), simplify((p(c)**2/(p(a)*p(b)**3)) / (c**2/(a*b**3))), simplify((p(c)**2/p(a)) / (c**2/a)))
assert check_item('THE-EQL-026')
