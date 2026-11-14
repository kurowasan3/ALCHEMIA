# THE-SOL-066 検証：混ざり合う2成分の有機層（理想溶液）と水の水蒸気蒸留
from sympy import Rational as Q, symbols, solve, N
P = 100  # 外圧 kPa
d90 = dict(W=Q(70), A=Q(56), B=Q(16))
d95 = dict(W=Q(845,10), A=Q(66), B=Q(195,10))
d125 = dict(A=Q(165), B=Q(65))
MA, MB, MW = 92, 106, 18
x = symbols('x')
xA = solve(d90['W'] + x*d90['A'] + (1-x)*d90['B'] - P, x)[0]
print('(1) x_A =', xA, N(xA))
assert xA == Q(7,20)
pA, pB, pW = xA*d90['A'], (1-xA)*d90['B'], d90['W']
print('(2) 分圧 W,A,B =', pW, N(pA), N(pB), ' 合計', N(pW+pA+pB))
yA_org = pA/(pA+pB)
print('    有機物中の A のモル分率', N(yA_org, 4))
mass_ratio = (pA*MA + pB*MB)/(pW*MW)
print('    有機物/水 質量比', N(mass_ratio, 5))
# (3) x_A = 0.20
for name, d in (('90', d90), ('95', d95)):
    tot = d['W'] + Q(1,5)*d['A'] + Q(4,5)*d['B']
    print('(3) x_A=0.20 の全蒸気圧', name, N(tot))
# 100 ℃ の水の蒸気圧 101 kPa > 100 なので必ず 100 ℃ 未満
# 単一成分 A のみ：90 ℃ で 70+56=126 > 100
print('    A のみ+水 90℃:', d90['W']+d90['A'])
# B のみ+水
print('    B のみ+水 90℃:', d90['W']+d90['B'], ' 95℃:', N(d95['W']+d95['B']))
# (4) 水なしで x_A=0.35 の有機混合物 125 ℃
tot125 = xA*d125['A'] + (1-xA)*d125['B']
print('(4) 125℃ 全圧', N(tot125))
assert tot125 == 100
yA_dry = xA*d125['A']/tot125
print('    水なし蒸気中の A のモル分率', N(yA_dry, 4))
print('    比 P_A/P_B 90℃', N(d90['A']/d90['B'],3), '125℃', N(d125['A']/d125['B'],3))
assert round(float(yA_org),2) == 0.65 and round(float(yA_dry),2) == 0.58
assert round(float(mass_ratio),1) == 2.3
print('OK')
