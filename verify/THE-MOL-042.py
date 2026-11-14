# THE-MOL-042 検証：物質量の単位を変えたときの定数・pH の変化（新単位 mol* = 粒子数が今の 1/10）
from sympy import Rational as Q, log, N, sqrt
f = Q(10)          # 同じ物質の量は mol* で 10 倍の数値になる（n* = 10 n）
def scale(power_mol):   # 単位に mol^k を含む量は 10^k 倍
    return f**power_mol
cases = {
 'Kc N2O4<=>2NO2 (mol/L)^1': scale(1),
 'Kc H2+I2<=>2HI (mol/L)^0': scale(0),
 'Kp N2O4<=>2NO2 Pa^1': scale(0),
 'Kw (mol/L)^2': scale(2),
 'Ksp Ag2CrO4 (mol/L)^3': scale(3),
 'F C/mol': scale(-1),
 'R Pa L/(K mol)': scale(-1),
 'Kf K kg/mol': scale(-1),
}
for k, v in cases.items(): print(k, v)
# 直接確認：N2O4 0.010 mol/L, NO2 0.020 mol/L
c1, c2 = Q('0.010'), Q('0.020')
assert (10*c2)**2/(10*c1) == 10*(c2**2/c1)
# (2) 中性の水
Kw = Q(1, 10**14); Kws = Kw*100
pH_neutral_star = -log(sqrt(Kws), 10)
print('pH* neutral =', N(pH_neutral_star)); assert abs(N(pH_neutral_star) - 6) < 1e-12
# pH* = 6.5 の溶液：[H+]* = 10^-6.5 mol*/L -> 現在の単位で 10^-7.5 mol/L
h_now = 10**Q('-6.5')/10
print('pH(now) =', N(-log(h_now, 10)), '-> 塩基性')
assert N(-log(h_now, 10)) > 7
# (3) 酢酸 0.10 mol/L, Ka = 2.7e-5
Ka = Q('2.7e-5'); c = Q('0.10')
a_now = sqrt(Ka/c); a_star = sqrt(10*Ka/(10*c))
print('alpha now', N(a_now), 'alpha*', N(a_star)); assert a_now == a_star
# 10 倍希釈の pH 差（強酸）：-log(c/10)-(-log c) = 1 はどちらの単位でも同じ
# 凝固点降下：Kf* m* = (Kf/10)(10 m)
Kf, m = Q('1.85'), Q('0.10'); assert (Kf/10)*(10*m) == Kf*m
# 変色域の pH：pH* = pH - 1
print('pH* = pH - 1 （例：フェノールフタレイン 8.0〜9.8 → 7.0〜8.8）')
print('OK')
