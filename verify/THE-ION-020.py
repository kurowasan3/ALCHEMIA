import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from math import log10
K1=4.0e-7; K2=5.0e-11
assert round(-log10(K1),1)==6.4 and round(-log10(K2),1)==10.3
h=4.0e-8; r1=K1/h; r2=K2/h; print(r1, r2, r1*r2); assert abs(r1-10)<1e-9 and abs(r1*r2-0.0125)<1e-12
h=10**-8.3; D=h*h+K1*h+K1*K2; print('frac HCO3 at 8.3', K1*h/D, h*h/D, K1*K2/D); assert K1*h/D>0.95
assert not check_item('THE-ION-020')
