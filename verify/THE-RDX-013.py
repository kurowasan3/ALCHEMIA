import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from sympy import Rational as R
so2=R(932,10000)/233*10; assert so2==R(4,1000)
i2rem=R(1,2)*R(1,10)*R(30,1000)*10; i2=R(5,100)*R(1,2); h2s=i2-i2rem-so2; print(so2,h2s); assert h2s==R(6,1000)
assert so2/h2s==R(2,3)
from chemlib import check_equation
assert check_equation('BaCl2 + H2SO4 -> BaSO4 + 2HCl') is None
assert not check_item('THE-RDX-013')
