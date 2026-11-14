import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, log, N
assert check_equation('4NH4+ + 6HCHO -> (CH2)6N4 + 4H+ + 6H2O') is None
c = R(1, 10)
nNH4 = c*R(146, 10)/1000*10          # 100 mL 中
N_pct = nNH4*14/1*100
AS_pct = nNH4/2*132/1*100
print('N%', N(N_pct, 6), 'AS%', N(AS_pct, 6))
assert round(float(N_pct), 1) == 20.4 and round(float(AS_pct), 1) == 96.4
Na = c*R(150, 10)/1000*10*14*100
print('(a) N%', N(Na, 5)); assert float(Na) > float(N_pct)
Ka = 5.7e-10
pKa = -log(Ka, 10)
frac = Ka/(Ka + 10**-8.5)
print('pKa', N(pKa, 4), 'pH8.5 で中和される割合', round(frac, 3))
H2SO4 = c*R(4, 10)/1000/2*10*98
print('遊離硫酸 g', N(H2SO4, 3))
print('OK')
