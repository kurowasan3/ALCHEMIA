import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from scipy.optimize import brentq
from math import sqrt, log10
K1=4.0e-7; K2=5.0e-11; Kw=1e-14
for C in (0.025, 2.5e-4):
    a=sqrt(K1/C); print('approx', C, a, -log10(C*a))
    def f(lh):
        h=10**lh; D=h*h+K1*h+K1*K2
        return h - Kw/h - C*(K1*h+2*K1*K2)/D
    h=10**brentq(f,-12,0); D=h*h+K1*h+K1*K2; co3=C*K1*K2/D
    print('exact pH', -log10(h), 'CO3', co3)
    assert round(-log10(h),1)==round(-log10(C*a),1) and abs(co3-5e-11)/5e-11<0.03
assert not check_item('THE-ION-013')
