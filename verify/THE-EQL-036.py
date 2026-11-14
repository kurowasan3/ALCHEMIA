# THE-EQL-036 検証：酢酸蒸気の二量体化（見かけの分子量と圧平衡定数）
from sympy import Rational as Q, symbols, solve, N, sqrt
R = Q(83,10)*1000; T = 400; V = Q(83,100); P = 10**5
n = P*V/(R*T); m = Q(240,100)
Mapp = m/n; print('n', n, 'M_app', Mapp); assert Mapp == 96
# M_app = 60 (pM + 2pD)/(pM + pD), pM + pD = P
pM, pD = symbols('pM pD', positive=True)
s = solve([pM + pD - P, 60*(pM + 2*pD) - Mapp*(pM + pD)], [pM, pD]); print(s)
Kp = s[pD]/s[pM]**2; print('Kp', Kp, N(Kp)); assert Kp == Q(3,80000)   # 3.75e-5 /Pa
# 全圧 2.0e4 Pa での見かけの分子量
p = symbols('p', positive=True)
pm = [r for r in solve(Kp*p**2 + p - 2*10**4, p) if r > 0][0]
pd = 2*10**4 - pm; M2 = 60*(pm + 2*pd)/(pm + pd); print('pM', pm, 'pD', pd, 'M', M2)
assert pm == Q(40000,3) and M2 == 80
# 見かけの分子量が 63 以下（二量体のモル分率 0.05 以下）となる全圧の上限
Pmax = Q(5,100)/(Kp*Q(95,100)**2); print('Pmax', N(Pmax)); assert abs(N(Pmax) - 1.48e3) < 5
# 二量体の物質量比（単量体単位の何 % が二量体か）at 1e5
frac = 2*s[pD]/(s[pM] + 2*s[pD]); print('二量体になっている割合', frac); assert frac == Q(3,4)
# 450 K で M=70 のとき Kp(450)
s3 = solve([pM + pD - P, 60*(pM + 2*pD) - 70*(pM + pD)], [pM, pD]); K450 = s3[pD]/s3[pM]**2
print('K450', N(K450)); assert K450 < Kp
print('OK')
