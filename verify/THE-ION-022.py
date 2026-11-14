import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from math import sqrt, log10
Kb=2e-5; Ksp=1.6e-11
oh=sqrt(Kb*0.10); Q=0.010*oh**2; print(oh, Q); assert Q>Ksp and abs(Q-2e-8)<1e-15
ohmax=sqrt(Ksp/0.010); nh4=Kb*0.10/ohmax; n=nh4*0.100; print(ohmax, nh4, n); assert abs(n-5e-3)<1e-12
pH=14+log10(ohmax); print(pH); assert round(pH,1)==9.6
assert 0.010*ohmax**3 > 1e-38; print(0.010*ohmax**3)
assert not check_item('THE-ION-022')
