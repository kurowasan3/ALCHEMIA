import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
import math
Ea = 2.3*8.3*0.90/(1/300-1/330); print(Ea); approx(Ea, 5.67e4, 0.005)
print('exact ln', 8.3*math.log(8)/(1/300-1/330))
e2 = 23e3/(2.3*8.3*300); e3 = 23e3/(2.3*8.3*330); print(e2, 10**e2, e3, 10**e3)
approx(10**e2, 1.0e4, 0.06); approx(10**e3, 4.5e3, 0.01)
assert check_item('THE-EQL-017')
