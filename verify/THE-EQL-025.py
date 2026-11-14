import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
from scipy.optimize import brentq
x = 1/3; approx(x/(1+x), 0.25, 1e-9)
approx(88*x, 29.3, 0.01)
V = (1+x)*8.3e3*500/1.0e5; approx(V, 55.3, 0.005)
Kc = (x/V)**2/((1-x)/V); approx(Kc, 3.0e-3, 0.01)
Kp = x**2/(1-x**2)*1.0e5
x2 = brentq(lambda t: t**2/(1-t**2)*2.0e5 - Kp, 1e-9, 1-1e-9); print('x at 2e5 Pa', x2); assert x2 < x
assert check_item('THE-EQL-025')
