import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R
for e in ['C6H12O6 -> 2C2H5OH + 2CO2','C6H12O6 -> 2C3H6O3','C2H5OH + 3O2 -> 2CO2 + 3H2O','C3H6O3 + 3O2 -> 3CO2 + 3H2O','C6H12O6 + 6O2 -> 6CO2 + 6H2O']:
    assert check_equation(e) is None, e
g,e,l=-2800,-1370,-1340
alc=g-2*e; lac=g-2*l; print(alc,lac); assert alc==-60 and lac==-120
a=R(32*31,2800)*100; b=R(2*31,120)*100; print(float(a),float(b)); assert round(a)==35 and round(b)==52
assert 14-31==-17
print(2800/120, 32/2)
print('THE-THC-030 OK')
