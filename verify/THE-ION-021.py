import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from math import log10, sqrt
S=lambda pH: 1e-22/10**(-2*pH)
KCu, KZn = 1e-30, 2e-18
s1=S(1.0); assert 0.01*s1>KCu and 0.01*s1<KZn; print('Cu left', KCu/s1)
h=sqrt(1e-22/(KZn/0.01)); pHmax=-log10(h); print('pHmax', pHmax); assert round(pHmax,2)==3.15
print('Cu at pHmax', KCu/(KZn/0.01)); assert abs(KCu/(KZn/0.01)-5e-15)<1e-20
assert 0.01*S(7.0)>KZn
assert not check_item('THE-ION-021')
