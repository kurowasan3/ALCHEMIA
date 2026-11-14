# THE-ION-103 検証：ギ酸の滴定から Ka と中和点の pH
from sympy import Rational, sqrt, log, N, nsolve, symbols
c = Rational(1, 10)*20/10  # 0.200 mol/L
H = 10**(-3.30)
print('c', c, '[H+]', H)
Ka = H*(0.5e-3)/(1.5e-3)  # [H+]×n_salt/n_acid
Ka_r = 5.0e-4/3            # log2=0.30 で [H+]=5.0e-4 とした値
print('Ka', Ka, Ka_r)
assert round(Ka_r*1e4, 1) == 1.7
cs = 2.00e-3/30.0e-3
OH = (cs*1e-14/Ka_r)**0.5
import math
pH = 14 + math.log10(OH)
print('salt conc', cs, 'OH', OH, 'pH', pH)
assert abs(pH - 8.30) < 0.01
# 厳密解（電荷収支）での 5.00 mL 点の pH を確認（近似の妥当性）
from sympy import solve
h = symbols('h', positive=True)
Ct = 2.0e-3/15e-3; Na = 0.5e-3/15e-3
hs = [x for x in solve(h*(Ka_r + h) + Na*(Ka_r + h) - Ct*Ka_r, h) if x > 0][0]
print('exact pH at 5 mL with Ka', -math.log10(float(hs)))
print('OK')
