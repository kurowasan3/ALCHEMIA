# THE-SOL-096 一定圧力で沸騰する溶液：蒸気量・残る溶媒・沸点上昇の連成
from sympy import symbols, nsolve, Rational as R, N
P, Rg, Kb, T0 = R(1013, 10)*1000, R(831, 100)*1000, R(52, 100), 373
n_solute = R(360, 100)/180; assert n_solute == R(2, 100)
w0 = 100
# (1) 近似 T=373 K で蒸発した水
V = 100
n1 = P*V/(Rg*T0); m1 = 18*n1
print('(1) n, m =', N(n1, 5), N(m1, 5)); assert round(float(m1), 1) == 58.8
# (2) 残りの水と沸点
wr = w0 - m1; mol = n_solute/(wr/1000); dT = Kb*mol
print('(2) water left', N(wr, 5), 'molality', N(mol, 4), 'dT', N(dT, 4), 'T', N(T0 + dT, 6))
assert round(float(wr), 1) == 41.2 and round(float(dT), 2) == 0.25
# 厳密に連立（蒸気量の T も沸点で）して比較
T = symbols('T')
f = lambda t: (t - T0) - Kb*n_solute/((w0 - 18*P*V/(Rg*t))/1000)
Tex = nsolve(f(T), T, 373.3)
mex = 18*P*V/(Rg*Tex)
print('exact T', Tex, 'evaporated', mex, 'diff g', N(m1 - mex, 3))
assert abs(float(m1 - mex) - 0.04) < 0.005
# (3) 温度比 a
a = R(573)/(T0 + dT); a0 = R(573, 373)
print('(3) a =', N(a, 6), 'pure', N(a0, 6)); assert round(float(a), 3) == 1.535 and a < a0
# (4) 沸点 374.0 K になるときの蒸気の体積
mol4 = R(1)/Kb; wr4 = n_solute/mol4*1000; ev4 = w0 - wr4
n4 = ev4/18; V4 = n4*Rg*374/P
print('(4) water left', N(wr4, 4), 'evap', N(ev4, 4), 'V', N(V4, 4))
assert round(float(V4)) == 153
print('all OK')
