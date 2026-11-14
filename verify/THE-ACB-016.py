import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from scipy.optimize import brentq
from math import log10
Ka=2.0e-5; Kw=1e-14
cH=0.100*8.0/20.0; cA=0.100*12.0/20.0; print(cH, cA); assert abs(cH-0.04)<1e-12 and abs(cA-0.06)<1e-12
def pH(v):
    Vt=20+v; Cl=cH*20/Vt; CA=cA*20/Vt; Na=0.1*v/Vt
    return -brentq(lambda lh: 10**lh+Na-Kw/10**lh-Cl-CA*Ka/(Ka+10**lh), -14, 1)
print(pH(0), pH(14)); assert round(pH(0),1)==1.4 and round(pH(14),1)==4.7
assert round(2-log10(4),1)==1.4 and round(-log10(2e-5),1)==4.7
assert not check_item('THE-ACB-016')
