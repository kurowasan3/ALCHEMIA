# INO-PPT-026 検証：沈殿の変化（Ag2CrO4 -> AgCl -> AgI）とモール法
import sys, json
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, sqrt
D = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))['results']
assert D['Ag+']['CrO4']['f'] == 'Ag2CrO4' and D['Ag+']['CrO4']['color'] == '暗赤色'
assert D['Ag+']['Cl']['f'] == 'AgCl' and D['Ag+']['Cl']['color'] == '白'
# 25 ℃ の代表的な溶解度積（教科書・便覧の値）
Ksp_AgCl, Ksp_Ag2CrO4, Ksp_AgI = 1.8e-10, 1.1e-12, 8.3e-17
K_conv1 = Ksp_Ag2CrO4 / Ksp_AgCl**2      # Ag2CrO4 + 2Cl- <=> 2AgCl + CrO4^2-
K_conv2 = Ksp_AgCl / Ksp_AgI             # AgCl + I- <=> AgI + Cl-
print('K(Ag2CrO4->AgCl)=%.2e  K(AgCl->AgI)=%.2e' % (K_conv1, K_conv2))
assert K_conv1 > 1e3 and K_conv2 > 1e3
# モル溶解度でも Ag2CrO4 > AgCl > AgI
s = {'Ag2CrO4': (Ksp_Ag2CrO4 / 4) ** (1 / 3), 'AgCl': Ksp_AgCl ** 0.5, 'AgI': Ksp_AgI ** 0.5}
print(s); assert s['Ag2CrO4'] > s['AgCl'] > s['AgI']
# モール法：[Cl-]=0.030, [CrO4^2-]≈1e-3 のとき，沈殿に必要な [Ag+]
need_AgCl = Ksp_AgCl / 0.030; need_Ag2CrO4 = (Ksp_Ag2CrO4 / 1e-3) ** 0.5
print('Ag+ for AgCl %.1e, for Ag2CrO4 %.1e' % (need_AgCl, need_Ag2CrO4)); assert need_AgCl < need_Ag2CrO4
# Cl- が 99.9 % 沈殿した時点 ([Cl-]=3e-5) でも AgCl が先
assert Ksp_AgCl / 3e-5 < need_Ag2CrO4
# 計算
c = R(50, 1000) * R(120, 10) / R(200, 10)
print('[Cl-] =', float(c)); assert c == R(3, 100)
for eq in ['2Ag+ + CrO4^2- -> Ag2CrO4', 'Ag2CrO4 + 2Cl- -> 2AgCl + CrO4^2-', 'AgCl + I- -> AgI + Cl-', 'AgCl + 2NH3 -> [Ag(NH3)2]+ + Cl-']:
    assert check_equation(eq) is None, eq
print('INO-PPT-026 OK')
