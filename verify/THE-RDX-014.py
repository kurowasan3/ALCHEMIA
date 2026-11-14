import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from sympy import Rational as R
# chain: 2Mn(OH)2 + O2 -> 2MnO(OH)2 ; MnO(OH)2 -> I2 ; I2 -> 2 S2O3
ratio = 2*1*2; assert ratio==4
from chemlib import check_equation
assert check_equation('Mn^2+ + 2OH- -> Mn(OH)2') is None
s=R(2,100)*R(6,1000); o2=s/4; mgL=o2*32*1000*R(1000,200); print(s,o2,mgL); assert mgL==R(48,10)
assert not check_item('THE-RDX-014')
