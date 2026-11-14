import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from sympy import Rational as R, N
e=R(4825,100)*1000/(R(965,10)*1000); assert e==R(1,2)
neg=96*e/2; pos=64*e/2; assert neg==24 and pos==16
acid=350-e*98; sol=1000-e*98+e*18; w=acid/sol*100; print(N(w)); assert round(float(w),1)==31.4
assert (e*98-e*18)==neg+pos
assert not check_item('THE-ELC-011')
