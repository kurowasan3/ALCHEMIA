import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from math import log10
E=lambda cu,zn: 1.10+0.059/2*log10(cu/zn)
assert E(2,1)>E(1,1) and E(1,0.5)>E(1,1); print(E(2,1), E(1,0.5))
assert not check_item('THE-ELC-012')
