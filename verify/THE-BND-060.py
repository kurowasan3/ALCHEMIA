# THE-BND-060 検証：ドライアイス（CO2 の分子結晶，面心立方格子）の密度と気体との比較
import sympy as sp
NA = sp.Rational(6, 1) * 10**23
a = sp.Rational(56, 100) * 10**-7  # cm
M = 44
n_cell = 8*sp.Rational(1, 8) + 6*sp.Rational(1, 2)
assert n_cell == 4
rho = n_cell*M/(NA*a**3)
print('density', float(rho)); assert round(float(rho), 1) == 1.7
mol_per_cm3 = n_cell/(NA*a**3)
print('mol in 1 cm3', float(mol_per_cm3))
V_gas = mol_per_cm3*sp.Rational(224, 10)*1000  # cm3 (0 ℃, 1 atm)
print('gas volume from 1 cm3', float(V_gas)); assert round(float(V_gas), -1) == 850
P = mol_per_cm3*sp.Rational(82, 1000)*323/sp.Rational(1, 1000)  # atm at 50 ℃, V = 1 cm3
print('ideal pressure at 323 K', float(P)); assert f"{float(P):.1e}" == '1.0e+03'
# 分子の体積：1 分子あたりの空間 a^3/4 と CO2 分子の大きさ（約 0.3–0.5 nm）が同程度 → 反発が支配的 → Z > 1
print('space per molecule nm^3', float((a*1e7)**3/4))
