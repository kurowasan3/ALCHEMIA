import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
from scipy.optimize import fsolve
import math
Kc = (0.60/3.0)**2/(0.30/3.0); approx(Kc, 0.40)
y = (2.1 - math.sqrt(2.1**2 - 4*0.5625))/2; print('y', y)
approx((0.75-y)/0.75, 0.58, 0.01)
R=8.3e3; T=330
P = (1.5-y)*R*T/6.0; approx(P, 5.4e5, 0.01)
# 独立検証：各室別々に解く。未知数 yA(N2O4 in A), yB, VA
NA, NB = 0.60, 0.15   # N2O4 換算
def f(v):
    yA, yB, VA = v; VB = 6.0-VA
    nA = yA + 2*(NA-yA); nB = yB + 2*(NB-yB)
    return [(2*(NA-yA)/VA)**2/(yA/VA)-Kc, (2*(NB-yB)/VB)**2/(yB/VB)-Kc, nA/VA - nB/VB]
yA, yB, VA = fsolve(f, [0.25, 0.06, 4.5]); print(yA, yB, VA, VA/(6-VA), (NA-yA)/NA, (NB-yB)/NB)
approx(VA/(6-VA), 4.0, 1e-4); approx((NA-yA)/NA, (0.75-y)/0.75, 1e-4)
assert check_item('THE-EQL-030')
