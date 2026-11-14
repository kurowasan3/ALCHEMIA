import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from math import log10, sqrt
from scipy.optimize import brentq
K1=4.0e-7; K2=5.0e-11; Kw=1e-14
print(K2/K1); assert abs(K2/K1-1.25e-4)<1e-12
p=-log10(sqrt(K1*K2)); print(p); assert round(p,2)==8.35
for C in (0.10, 0.010):
    def f(lh):
        h=10**lh; D=h*h+K1*h+K1*K2
        return h + C - Kw/h - C*(K1*h+2*K1*K2)/D
    ph=-brentq(f,-14,0); print(C, ph); assert abs(ph-8.35)<0.05
assert not check_item('THE-ION-018')
