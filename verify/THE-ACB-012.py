import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
import math
from scipy.optimize import brentq
h = (0.1+0.001)/2; print(-math.log10(h)); approx(round(-math.log10(h),1), 1.3, 1e-9)
x = brentq(lambda v: 0.10*(v-20)/(20+v) - 0.010, 20, 100); print(x); approx(round(x,1), 24.4, 1e-9)
Ka = 2.5e-5; c = (1e-3)**2/Ka + 1e-3; c2 = c/10; y = (-Ka+math.sqrt(Ka*Ka+4*Ka*c2))/2; pH = -math.log10(y); print('acetic diluted pH', pH); assert 3.0 < pH < 4.0
assert check_item('THE-ACB-012')
