# THE-SOL-076 検証：フェニル酢酸の凝固点降下（ベンゼン中の二量体化と酢酸中で会合が見えない理由）
import sympy as sp
R = sp.Rational
M = 136
# ベンゼン：1.36 g / 50.0 g, Kf 5.12, ΔT 0.563
c0 = R(136, 100)/M/R(50, 1000)          # 0.200 mol/kg（すべて単量体なら）
assert c0 == R(1, 5)
mtot = R(563, 1000)/R(512, 100)
print('particles', float(mtot))
m, d = sp.symbols('m d', positive=True)
s = sp.solve([m + 2*d - c0, m + d - mtot], [m, d]); print(s)
frac = 2*s[d]/c0; K = s[d]/s[m]**2
print('dimer fraction of acid', float(frac), ' K', float(K), ' apparent M', float(M*c0/mtot))
assert round(float(frac), 2) == 0.90 and round(float(K), -1) == 230 or round(float(K)) == 225
# 酢酸溶媒：0.680 g / 50.0 g, Kf 3.90, 観測 ΔT 0.390 → 見かけの分子量
Mapp = R(680, 1000)/(R(390, 1000)/R(390, 100))/R(50, 1000)
print('apparent M in acetic acid', Mapp); assert Mapp == 136
