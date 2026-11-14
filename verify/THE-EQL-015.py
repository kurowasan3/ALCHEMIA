import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
from sympy import symbols, solve, Rational, log
k, C0, t = symbols('k C0 t', positive=True)
def tt(frac): return solve(1/(frac*C0) - 1/C0 - k*t, t)[0]
print(tt(Rational(1,2)), tt(Rational(1,4)))
kc = Rational(1, 50)  # kC0
approx(float(tt(Rational(1,4)).subs(k, kc/C0)), 150)
approx(2*50, 100)
x = symbols('x', positive=True); kc2 = solve(1/Rational(3,4) - 1 - x*60, x)[0]; print('k\'C0', kc2)
approx(float(1/kc2), 180)
assert check_item('THE-EQL-015')
