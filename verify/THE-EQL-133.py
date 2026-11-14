import sys, sympy as sp
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
K = 4.0  # 1e5 Pa ; 0.100 mol <-> 1.00e5 Pa at V0
y = sp.symbols('y', positive=True)
y0 = [s for s in sp.solve((2*y)**2 - K*(2 - y), y) if s < 2][0]
assert abs(y0 - 1) < 1e-9; print('min C g', float(y0*0.1*12))
nC = 0.72/12; nCO2 = 0.200 - nC; nCO = 2*nC
Q = (nCO/0.1)**2/(nCO2/0.1); print('Q', Q); assert Q < K
k = K/Q; print('k', k, 'V/V0', 1/k); assert round(1/k, 2) == 0.26
f = 8; c2 = nCO2/0.1*f; c = nCO/0.1*f
z = sp.symbols('z')
zs = [s for s in sp.solve((c - 2*z)**2 - K*(c2 + z), z) if 0 < s < c/2]
assert len(zs) == 1
g = float(zs[0])/f*0.1*12; print('z', float(zs[0]), 'C g', g); assert round(g, 2) == 0.19
for w in [0, 0.6, 1.2]: print(w, 2.00 + w/1.2)
assert check_equation('C + CO2 <=> 2CO') is None and check_equation('2CO -> C + CO2') is None
print('OK')
