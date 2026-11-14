# INO-TRN-037 検証：亜硫酸水素イオンによる六価クロムの還元と水酸化クロム(III)の沈殿
import sys, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
eqs = ['Cr2O7^2- + 14H+ + 6e- -> 2Cr^3+ + 7H2O', 'HSO3- + H2O -> SO4^2- + 3H+ + 2e-',
       'Cr2O7^2- + 3HSO3- + 5H+ -> 2Cr^3+ + 3SO4^2- + 4H2O', 'Cr^3+ + 3OH- -> Cr(OH)3']
for e in eqs: r = check_equation(e); print(e, r); assert r is None
Cr = 52.0e-3/52            # mol（排水 1.00 L 中の Cr）
nCr2O7 = Cr/2; nHSO3 = 3*nCr2O7; mNaHSO3 = nHSO3*104*1e3
nH = 5*nCr2O7
print('Cr2O7', nCr2O7, ' NaHSO3', nHSO3, 'mol', mNaHSO3, 'mg', ' 消費 H+', nH)
assert abs(mNaHSO3 - 156) < 1e-9
Ksp = 6.3e-31
for pH in (5.0, 8.0):
    OH = 10**(pH - 14); c = Ksp/OH**3; print('pH', pH, '[Cr3+] =', c, 'mol/L =', c*52e3, 'mg/L')
# 基準 0.52 mg/L = 1.0e-5 mol/L となる pH
OH = (Ksp/1.0e-5)**(1/3); pH = 14 + math.log10(OH); print('最低 pH', pH)
assert abs(round(pH, 1) - 5.6) < 1e-9
print('log10 6.3 =', math.log10(6.3))
print('OK')
