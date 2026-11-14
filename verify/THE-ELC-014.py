import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from sympy import Rational as R, N
w=R(36,1)/18; e=2*w; Q=e*96500; W=Q*R(9,10); print(Q, N(W/1000)); assert round(float(W/1000))==347
eff=W/1000/(286*w)*100; print(N(eff)); assert round(float(eff))==61
assert not check_item('THE-ELC-014')
