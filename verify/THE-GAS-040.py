# THE-GAS-040 検証：ピストンで仕切った 2 室の一方だけを加熱
from sympy import symbols, solve, Rational as Q, N
P0 = Q(10)**5; L0 = 30; T0 = 300; T1 = 360   # 27 ℃ → 87 ℃
# (1) B 室が大気に通じている：A は定圧
LA1 = L0*Q(T1, T0); print('(1) LA =', LA1, ' 移動', LA1-L0, 'cm')
assert LA1 - L0 == 6
# (2) 両方閉じている
LA, P = symbols('LA P', positive=True)
sol = solve([P*LA/T1 - P0*L0/T0, P*(60-LA)/T0 - P0*L0/T0], [LA, P], dict=True)[0]
print('(2) LA =', N(sol[LA],5), ' 移動', N(sol[LA]-L0,3), ' P =', N(sol[P],4))
assert abs(float(sol[LA]-L0) - 2.7) < 0.05 and abs(float(sol[P]) - 1.1e5) < 0.05e5
# B 室の圧力も A と同じ (P) ・ B 室の温度は不変
# (4) 両方閉じたまま A を何 K にすると 6.0 cm 動くか
PB = P0*L0/24; TA = T0*PB*36/(P0*L0)
print('(4) P =', N(PB,4), ' TA =', TA, 'K =', TA-273, '℃')
assert TA == 450
print('OK')
