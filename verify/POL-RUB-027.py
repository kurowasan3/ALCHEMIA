# POL-RUB-027 検証：陽イオン交換樹脂のバッチ平衡
from sympy import symbols, solve
cap = 4.00e-3*10.0; nNa0 = 0.100*0.500
H_rel = 0.100*30.0e-3/50.0e-3*0.500
RH = cap - H_rel; RNa = H_rel; Na = nNa0 - H_rel
K = RNa*H_rel/(RH*Na)
print('H released', H_rel, 'RH', RH, 'Na sol', Na, 'K', K)
assert abs(H_rel-0.030) < 1e-12 and abs(K-4.5) < 1e-9
x = symbols('x', positive=True)
sols = [s for s in solve(x**2 - K*(0.080-x)*(0.050-x), x) if 0 < s < 0.050]
frac = float(sols[0])/0.050; print('x', sols, 'frac', frac)
assert round(frac*100) == 81
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
assert check_equation('NaOH + HCl -> NaCl + H2O') is None
print('OK')
