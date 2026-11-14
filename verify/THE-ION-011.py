import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
import math
from scipy.optimize import brentq
c = 0.100*16.0/10.0; approx(c, 0.160)
Ka = 10**-4.70; approx(Ka, 2.0e-5, 0.003)
approx(4.70 - 0.48, 4.22, 1e-9)
def pH(vb, water, Ka=2.0e-5):
    V = 10.0+water+vb; Ca = c*10.0/V; Na = 0.100*vb/V
    f = lambda p: (Na + 10**-p - 1e-14/10**-p) - Ca*Ka/(Ka+10**-p)
    return brentq(f, 0, 14)
for w in (40.0, 0.0):
    print('water', w, 'pH0', pH(0, w), 'pH4', pH(4.0, w), 'pH8', pH(8.0, w))
assert abs(pH(8.0, 40.0) - 4.70) < 0.01 and abs(pH(8.0, 0.0) - 4.70) < 0.01
assert abs(pH(4.0, 40.0) - 4.22) < 0.01
assert pH(0, 0.0) < pH(0, 40.0)
assert check_item('THE-ION-011')
