"""INO-HAL-014 検証：海藻灰からのヨウ素の取り出し（反応式のつり合いと量的関係）"""
import sys
from sympy import Rational as R
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for eq in ['2NaI + MnO2 + 2H2SO4 -> Na2SO4 + MnSO4 + I2 + 2H2O', '2I- + MnO2 + 4H+ -> Mn^2+ + I2 + 2H2O',
           '2I- + H2O2 + 2H+ -> I2 + 2H2O', 'MnO2 + 4H+ + 2e- -> Mn^2+ + 2H2O', '2I- -> I2 + 2e-',
           'H2O2 + 2H+ + 2e- -> 2H2O', 'Cl2 + 2I- -> 2Cl- + I2']:
    r = check_equation(eq); print(eq, r or 'OK'); assert r is None
nI = R(4, 100)                     # NaI 0.0400 mol
MnO2 = 55 + 2 * 16; I2 = 2 * 127; H2O2 = 2 * R(1) + 2 * 16
mMn = nI / 2 * MnO2; mI2 = nI / 2 * I2; mH2O2sol = nI / 2 * H2O2 / R(3, 100)
print('MnO2 g', float(mMn), 'I2 g', float(mI2), '3.0% H2O2 水 g', float(mH2O2sol))
assert round(float(mMn), 2) == 1.74 and round(float(mI2), 2) == 5.08 and round(float(mH2O2sol), 1) == 22.7
# MnO2・H2O2 はどちらも 1 mol で電子 2 mol を受け取るので必要量は I2 と同じ物質量
print('INO-HAL-014 all OK')
