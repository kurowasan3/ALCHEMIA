import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from sympy import Rational as R
nox=R(670,1000)/134*R(25,250); assert nox==R(5,10000)
c=R(2,5)*nox/R(20,1000); print(c); assert c==R(1,100)
nh=R(5,2)*c*R(176,10000); print(nh); assert nh==R(44,100000)
conc=nh/(R(5,1000)*R(10,100)); print(float(conc)); assert conc==R(88,100)
print('mass% approx', float(conc*34/1010*100))
assert not check_item('THE-RDX-011')
