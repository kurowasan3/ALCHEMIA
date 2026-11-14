# INO-AMP-013 検証：SnCl2 と FeCl3 の量的関係・反応式・検出反応
import sys, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R
for eq in ['Sn^2+ + 2Fe^3+ -> Sn^4+ + 2Fe^2+', 'Sn + Sn^4+ -> 2Sn^2+', '2Sn^2+ + O2 + 4H+ -> 2Sn^4+ + 2H2O']:
    assert check_equation(eq) is None, eq
Fe = R(20, 1000) * R(50, 1000)              # mol
need_Sn = Fe / 2
V_need = need_Sn / R(10, 1000) * 1000       # mL
assert V_need == 50, V_need
for V, label in [(30, 'X'), (80, 'Y')]:
    Sn = R(10, 1000) * R(V, 1000)
    reduced = min(2 * Sn, Fe)
    Fe3_left = Fe - reduced
    Sn_left = Sn - reduced / 2
    print(label, 'Fe3+ left =', float(Fe3_left), 'Sn2+ left =', float(Sn_left))
    if label == 'X': assert Fe3_left == R(4, 10000)
    if label == 'Y': assert Fe3_left == 0 and Sn_left == R(3, 10000)
ions = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))
print('Fe3+ SCN:', ions['results']['Fe^3+']['SCN'])
print('Fe2+ SCN:', ions['results']['Fe^2+']['SCN'])
assert ions['results']['Fe^2+']['SCN']['kind'] == 'none'
print('INO-AMP-013 OK')
