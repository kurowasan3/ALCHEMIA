# INO-HAL-021 検証：銅法・ブローアウト法の反応式と量的関係
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, nsimplify
eqs = ['2Cu^2+ + 4I- -> 2CuI + I2',
       'Cu^2+ + Fe^2+ + I- -> CuI + Fe^3+',
       '2CuI + O2 -> 2CuO + I2',
       'Cl2 + 2I- -> 2Cl- + I2',
       'I2 + 5Cl2 + 6H2O -> 2IO3- + 10Cl- + 12H+',
       'I2 + SO2 + 2H2O -> 2I- + SO4^2- + 4H+',
       'ClO- + 2I- + 2H+ -> Cl- + I2 + H2O',
       '2NaI + 2CuSO4 + 2FeSO4 -> 2CuI + Fe2(SO4)3 + Na2SO4']
for e in eqs:
    r = check_equation(e); print(e, '->', r); assert r is None
# 地下水 5.0 m3, I- 127 mg/L
nI = R(5000)*R(127,1000)/127
print('n(I-)=', nI); assert nI == 5
# (1) CuSO4 だけ：沈殿に入る I の割合 1/2
frac = R(2,4); print('割合', frac)
# (2) 銅法：CuSO4, FeSO4 各 5.00 mol, O2 2.50 mol, I2 2.50 mol = 635 g
nI2 = nI/2; print('I2', nI2, 'mol', float(nI2*254), 'g'); assert nI2*254 == 635
print('O2', nI2)
# (3) Cl2 4.00 mol 吹き込み：2.50 mol で I2 2.50 mol、残り 1.50 mol が I2 を IO3- に（I2:Cl2=1:5）
excess = R(4) - nI/2
lost = excess/5
left = nI2 - lost
print('I2 残り', left, float(left)); assert left == R(22,10)
# (4) SO2・NaClO：I2 2.20 mol -> SO2 2.20 mol、NaClO 2.20 mol
print('SO2 =', left, 'NaClO =', left)
# 酸化数：S +4 -> +6, I 0 -> -1, Cl(ClO-) +1 -> -1
# (5) ブローアウト法で I2 1 mol あたり Cl2 1, SO2 1, NaClO 1（理想）
# 電子の数の確認：Cl2 1 mol = 2 e-, I2->2I- 2e-, ClO- 2 e-
print('OK')
