# THE-RDX-045 検証：液性で変わる MnO4- の還元生成物を使った Fe2+・Mn2+ の同時定量
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, Eq
c = 0.0200; V1 = 12.0e-3; V2 = 30.0e-3
x, y = symbols('x y', positive=True)
sol = solve([Eq(5*c*V1, y), Eq(3*c*V2, y + 2*x)], [x, y])
print(sol)
X, Y = float(sol[x]), float(sol[y])
assert abs(Y-1.20e-3) < 1e-15 and abs(X-3.0e-4) < 1e-15
print('Fe2+', Y/0.0100, 'Mn2+', X/0.0100)
mMnO2 = (X + c*V2)*(54.9+32.0)*1000; print('MnO2 mg', mMnO2); assert round(mMnO2) == 78
xw = (5*c*V2 - Y)/2; print('wrong Mn2+', xw/0.0100); assert abs(xw/0.0100 - 0.0900) < 1e-12
for e in ['MnO4- + 8H+ + 5e- -> Mn^2+ + 4H2O', 'MnO4- + 4H+ + 3e- -> MnO2 + 2H2O', 'Mn^2+ + 2H2O -> MnO2 + 4H+ + 2e-',
          '2MnO4- + 3Mn^2+ + 2H2O -> 5MnO2 + 4H+', 'MnO4- + 5Fe^2+ + 8H+ -> Mn^2+ + 5Fe^3+ + 4H2O',
          'MnO4- + 3Fe^2+ + 4H+ -> MnO2 + 3Fe^3+ + 2H2O', 'ZnO + 2H+ -> Zn^2+ + H2O', 'Fe^2+ -> Fe^3+ + e-']:
    assert check_equation(e) is None, e
print('OK')
