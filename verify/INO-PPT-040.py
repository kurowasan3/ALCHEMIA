# INO-PPT-040 検証：アンモニア水・水酸化ナトリウム水溶液による Al・Fe・Zn の分離と重量分析
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
import sympy as sp
R = sp.Rational
MFe2O3, MAl2O3 = 2*56 + 3*16, 2*27 + 3*16
assert (MFe2O3, MAl2O3) == (160, 102)
mC, mD = R(364, 1000), R(160, 1000)
nFe = 2*mD/MFe2O3; nAl = 2*(mC - mD)/MAl2O3
print('Fe', nFe, 'Al', nAl); assert nFe == R(2, 1000) and nAl == R(4, 1000)
cFe, cAl = nFe/R(1, 10), nAl/R(1, 10)
print('conc', float(cFe), float(cAl))
for e in ['Al^3+ + 3NH3 + 3H2O -> Al(OH)3 + 3NH4+', 'Fe^3+ + 3NH3 + 3H2O -> Fe(OH)3 + 3NH4+',
          'Zn^2+ + 4NH3 -> [Zn(NH3)4]^2+', 'Al^3+ + 4OH- -> [Al(OH)4]-', 'Zn^2+ + 4OH- -> [Zn(OH)4]^2-',
          'Fe^3+ + 3OH- -> Fe(OH)3', '2Fe(OH)3 -> Fe2O3 + 3H2O', '2Al(OH)3 -> Al2O3 + 3H2O']:
    assert check_equation(e) is None, e
print('equations OK')
