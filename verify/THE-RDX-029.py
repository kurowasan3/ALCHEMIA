# THE-RDX-029 検証：亜硫酸ナトリウムを含む硫酸鉄(II)七水和物の KMnO4 滴定（煮沸の有無）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from sympy import symbols, solve
from chemlib import check_equation
H, O, Na, S, Fe = 1.0, 16, 23, 32, 56
M_F = Fe + S + 4*O + 7*(2*H + O)
M_S = 2*Na + S + 3*O
print('M', M_F, M_S)
c, VA, VB = 0.0200, 20.0e-3, 24.0e-3
f, s = symbols('f s')
sol = solve([f - 5*c*VA, f + 2*s - 5*c*VB], [f, s])
print(sol)
nf, ns = float(sol[f]), float(sol[s])
assert abs(nf - 2.00e-3) < 1e-12 and abs(ns - 2.00e-4) < 1e-12
W = 5*nf*M_F + 5*ns*M_S
print('sample mass', W)
assert abs(W - 2.906) < 1e-9
pur = 5*nf*M_F/W*100
app = 5*(5*c*VB)*M_F/W*100
print('purity', pur, 'apparent', app)
assert round(pur, 1) == 95.7 and round(app) == 115
for eq in ['SO3^2- + 2H+ -> SO2 + H2O', 'MnO4- + 5Fe^2+ + 8H+ -> Mn^2+ + 5Fe^3+ + 4H2O',
           '2MnO4- + 5SO2 + 2H2O -> 2Mn^2+ + 5SO4^2- + 4H+', '4Fe^2+ + O2 + 4H+ -> 4Fe^3+ + 2H2O']:
    r = check_equation(eq); print(eq, r); assert r is None
print('OK')
