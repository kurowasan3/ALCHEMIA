import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
import math
kf1 = 4.0e-3/(0.5*0.4); kf2 = 1.2e-2/(0.5*0.4); kr1 = 8.0e-5/0.04; kr2 = 4.8e-4/0.04
for a, b in [(kf1, 2.0e-2), (kr1, 2.0e-3), (kf2, 6.0e-2), (kr2, 1.2e-2)]: approx(a, b, 1e-9)
approx(kf1/kr1, 10, 1e-9); approx(kf2/kr2, 5.0, 1e-9)
# アレニウス：倍率 ratio = exp(Ea/R (1/T1-1/T2)) → 倍率大ほど Ea 大
x = (1/300 - 1/320)
Eaf = math.log(3)*8.314/x; Ear = math.log(6)*8.314/x; print('Ea f, r (例 T1=300,T2=320)', Eaf, Ear); assert Ear > Eaf
assert check_item('THE-EQL-028')
