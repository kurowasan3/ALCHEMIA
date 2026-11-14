import sys, sympy as sp
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
P1, P2 = 1.00, 1.60   # 1e5 Pa
x = P2 - P1; pCO2 = P1 - x; pCO = 2*x; K = pCO**2/pCO2
print(pCO2, pCO, K); assert abs(K-3.6) < 1e-9
y = sp.symbols('y', positive=True)
sol = [s for s in sp.solve((2*pCO - 2*y)**2 - K*(2*pCO2 + y), y) if 0 < s < pCO]
assert len(sol) == 1; y0 = float(sol[0])
tot = 2*pCO2 + y0 + 2*pCO - 2*y0
print('y', y0, 'total', tot); assert round(tot, 1) == 3.0 and tot < 3.2
assert check_equation('C + CO2 <=> 2CO') is None
print('OK')
