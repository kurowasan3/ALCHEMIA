# THE-THC-059 溶解熱と溶解度の温度変化が連成する断熱溶解（場合分け・しきい値・水の追加）
from sympy import symbols, solve, Rational as R, N
T, W = symbols('T W')
dH, M = R(24000), R(100)            # J/mol, g/mol
s = lambda t: 2*t - 10              # g/水100 g（0〜40 ℃）
C0 = R(480)                         # 水 100 g を含む容器と内容物の熱容量 J/K
q_per_g = dH/M                      # 240 J/g
T0 = 30
# (1) 20.0 g
Tf = T0 - q_per_g*20/C0
print('(1) T =', Tf, ' s =', s(Tf)); assert Tf == 20 and s(Tf) >= 20
# (2) 40.0 g: 全部溶けると仮定
Tall = T0 - q_per_g*40/C0
print('(2) all-dissolve T =', Tall, 's =', s(Tall)); assert s(Tall) < 40
T2 = solve(C0*(T0 - T) - q_per_g*s(T), T)[0]
m2 = s(T2)
print('(2) saturated T =', T2, 'm =', m2); assert T2 == R(35, 2) and m2 == 25
# (3) しきい値 w*: 全部溶けて丁度飽和
w = symbols('w')
ws = solve(s(T0 - q_per_g*w/C0) - w, w)[0]
print('(3) w* =', ws); assert ws == 25
# (4) 状態 (2) に 17.5 ℃ の水 100 g（比熱 4.0）を加える
C1 = C0 + 400
T4 = solve(C1*(T2 - T) - q_per_g*(2*s(T) - m2), T)[0]
m4 = 2*s(T4)
print('(4) T =', N(T4, 5), ' dissolved total =', N(m4, 5)); assert m4 < 40 and m4 > m2
assert round(float(T4), 1) == 14.2 and round(float(m4), 1) == 37.0
# 全部溶けると仮定した場合の矛盾も確認
Tall4 = T2 - q_per_g*(40 - m2)/C1
assert 2*s(Tall4) < 40
# (5) 40.0 g をすべて溶かすのに必要な 17.5 ℃ の水の追加量の最小値
eqs = [ (100 + W)/100*s(T) - 40, (C0 + 4*W)*(T2 - T) - q_per_g*(40 - m2) ]
sol = [x for x in solve(eqs, [W, T], dict=True) if x[W] > 0]
print('(5)', [(N(x[W], 5), N(x[T], 5)) for x in sol])
print('all OK')
