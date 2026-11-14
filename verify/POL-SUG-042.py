# POL-SUG-042 検証：デンプンの加水分解とアルコール発酵（質量減少からの消費率）
import sys
from sympy import Rational as Q, N
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for eq in ['C6H12O6 -> 2C2H5OH + 2CO2', 'C12H22O11 + H2O -> 2C6H12O6',
           'C6H10O5 + H2O -> C6H12O6']:
    assert check_equation(eq) is None, eq
H, C, O, Cu = Q(1), Q(12), Q(16), Q('63.5')
M_unit = 6*C + 10*H + 5*O      # 162
M_glc = 6*C + 12*H + 6*O       # 180
M_mal = 12*C + 22*H + 11*O     # 342
M_co2 = C + 2*O                # 44
M_et = 2*C + 6*H + O           # 46
M_cu2o = 2*Cu + O
starch = Q(81)
n_glc = starch / M_unit
print('glucose mol', n_glc, 'mass', n_glc*M_glc)
assert n_glc == Q(1,2) and n_glc*M_glc == 90
loss = Q('17.6')
n_co2 = loss / M_co2
n_used = n_co2 / 2
pct = n_used / n_glc * 100
print('CO2', n_co2, 'used', n_used, 'pct', pct)
assert pct == 40
et = n_co2 * M_et
print('ethanol', et); assert et == Q('18.4')
rest = n_glc - n_used
cu2o = rest * M_cu2o
print('Cu2O', N(cu2o, 4)); assert round(float(cu2o), 1) == 42.9
n_from_mal = 2 * starch / M_mal
print('maltose -> glucose mol', N(n_from_mal, 4))
assert n_from_mal < n_glc and round(float(n_from_mal), 3) == 0.474
print('OK')
