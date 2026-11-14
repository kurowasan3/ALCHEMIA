import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
from sympy import Rational as R
data = [(0.020, 1.2e-4), (0.040, 2.4e-4), (0.060, 3.6e-4)]
ks = [v/c for c, v in data]; print('k', ks)
for k in ks: approx(k, 6.0e-3, 1e-9)
v = 6.0e-3*0.050
approx(v*6/4, 4.5e-4, 1e-9); approx(v/4, 7.5e-5, 1e-9)
assert check_item('THE-EQL-014')
