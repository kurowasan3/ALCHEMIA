# POL-RES-035 検証：結合の向きに相関がある鎖（自由回転鎖）の末端間距離
from sympy import symbols, summation, simplify, Rational as R, sqrt, log, nsimplify, solve, cos, acos, pi, N as Nn
a, n, l, Nb = symbols('a n l N', positive=True)
# (1) 公式の確認：<R^2>/b^2 = N + 2 sum_{n=1}^{N-1} sum_{l=1}^{N-n} a^l
def R2(av, Nv):
    s = Nv
    for i in range(1, Nv):
        for j in range(1, Nv - i + 1):
            s += 2 * av**j
    return s
def formula(av, Nv):
    return Nv*(1+av)/(1-av) - 2*av*(1-av**Nv)/(1-av)**2
for av in [R(1,3), R(98,100), R(4,5)]:
    for Nv in [1, 2, 3, 5, 8, 20]:
        assert simplify(R2(av, Nv) - formula(av, Nv)) == 0
# 記号でも確認
F = Nb*(1+a)/(1-a) - 2*a*(1-a**Nb)/(1-a)**2
print('公式 OK')
# (2) a = 1/3（四面体の結合角 109.5° -> 結合ベクトルのなす角 70.5°，cos = 1/3）
print('cos(180-109.47)=', float(cos(pi - acos(-R(1,3)))))
a1 = R(1,3); Cinf = (1+a1)/(1-a1); print('C_inf', Cinf); assert Cinf == 2
CN = lambda av, Nv: formula(av, Nv)/Nv
for Nv in range(1, 12):
    print(Nv, float(CN(a1, Nv)), CN(a1, Nv) >= R(9,10)*Cinf)
Nmin = min(Nv for Nv in range(1, 100) if CN(a1, Nv) >= R(9,10)*Cinf); print('Nmin', Nmin); assert Nmin == 8
print('C_10', float(CN(a1, 10)))
# (3) a = 0.98：N = 5 と N = 5000 の sqrt(<R^2>)/b，棒（N）・C_inf N の平方根との比較
a2 = R(98,100); Cinf2 = (1+a2)/(1-a2); print('C_inf2', Cinf2); assert Cinf2 == 99
import math
ff = lambda av, Nv: Nv*(1+av)/(1-av) - 2*av*(1-av**Nv)/(1-av)**2   # 浮動小数点版
for Nv in [5, 5000]:
    r = math.sqrt(ff(0.98, Nv)); print(Nv, r, 'rod', Nv, 'coil', math.sqrt(99*Nv))
r5 = math.sqrt(ff(0.98, 5)); r5000 = math.sqrt(ff(0.98, 5000))
# 傾き：N=5 付近（5 -> 10）と N=5000 付近（5000 -> 10000）
for (N1, N2) in [(5, 10), (5000, 10000)]:
    s_ = (math.log10(math.sqrt(ff(0.98, N2))) - math.log10(math.sqrt(ff(0.98, N1)))) / math.log10(N2/N1)
    print('傾き', N1, N2, s_)
# (4) 測定：b = 0.150 nm，N = 2.00e4，sqrt<R^2> = 63.6 nm -> C = 9.0 -> a = 0.80
C = (R('63.6')**2) / (R('2.00e4') * R('0.150')**2); print('C', float(C))
av = solve((1+a)/(1-a) - R(9), a)[0]; print('a', av); assert av == R(4,5)
# 有限長の補正の大きさ
corr = 2*float(av)*(1-0.8**20000)/(20000*(1-0.8)**2); print('補正項', float(corr), '相対', corr/9)
print('角度', float(acos(R(4,5))*180/pi))
# (5) ポリエチレンの実測 C ≈ 7 は自由回転の 2 より大きい
print('OK')
