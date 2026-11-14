import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from sympy import Rational as R, N
c1=R(32,1000)/32*R(3,2)*32*1000; assert c1==48
c2=R(32,1000)/60*2*32*1000; print(N(c2)); assert round(float(c2),1)==34.1
e=5*R(5,1000)*R(6,1000); o2=e/4; cod=o2*32*1000*10; print(e, o2, cod); assert cod==12
assert not check_item('THE-RDX-012')
