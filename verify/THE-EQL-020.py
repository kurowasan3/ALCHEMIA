import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
from sympy import symbols, solve
R=8.3e3; T=350; n0=0.0200
n = 9.30e4*1.00/(R*T); a = n/n0-1; print('alpha', a); approx(round(a,2), 0.60, 1e-9)
Kp = 4*0.6**2/(1-0.6**2)*9.30e4; approx(Kp, 2.09e5, 0.005)
Kc = 4*n0*0.6**2/(0.4*1.00); approx(Kc, 0.072, 1e-6)
x = symbols('x', positive=True)
a2 = [s for s in solve(4*n0*x**2/((1-x)*0.50) - Kc, x) if s < 1][0]; print('alpha2', a2); approx(float(a2), 0.48, 0.01)
P2 = n0*(1+float(a2))*R*T/0.50; print(P2); approx(P2, 1.7e5, 0.02); assert P2 < 2*9.30e4
assert check_item('THE-EQL-020')
