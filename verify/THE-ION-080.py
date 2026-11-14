# THE-ION-080 検証：蒸発濃縮による硫酸カルシウム二水和物の析出
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, sqrt, Rational as Q
assert check_equation('Ca^2+ + SO4^2- + 2H2O -> CaSO4.2H2O') is None
assert check_equation('CaSO4.2H2O -> CaSO4 + 2H2O') is None
Ca0, S0, K = Q(15, 1000), Q(40, 1000), Q(24, 10000)   # mol/kg, (mol/kg)^2
X = symbols('X', positive=True)
Xs = solve(X**2*Ca0*S0 - K, X)[0]; print('(1) X =', Xs); assert Xs == 2
# (2) 水を 1/4 (0.25 kg) に濃縮：析出量 s mol（もとの水 1 kg あたり）
s = symbols('s')
sol = [r for r in solve((Ca0 - s)*(S0 - s) - K*Q(1, 4)**2, s) if 0 < r < Ca0]
print('(2) s =', sol); assert sol == [Q(1, 100)]
M = 40.1 + 32.1 + 4*16.0 + 2*18.0
print('   mass =', M*0.010, ' fraction of Ca =', float(sol[0]/Ca0))
assert round(M*0.010, 2) == 1.72
# (3) X=2 と X=4 での濃度
for x, ss in [(2, 0), (4, Q(1, 100))]:
    mc, ms = (Ca0 - ss)*x, (S0 - ss)*x
    print('X =', x, ' Ca =', float(mc), ' SO4 =', float(ms), ' product =', float(mc*ms))
    assert mc*ms == K
# 析出後も (SO4 - Ca) の物質量 0.025 mol は保存される
print('OK')
