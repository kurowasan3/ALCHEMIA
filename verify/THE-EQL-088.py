# THE-EQL-088 検証：金属酸化物を使う2段階の水の熱化学分解サイクル（log Kp = C1/T + C2）
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from sympy import Rational as Q, symbols, solve, N, log
from chemlib import check_equation
T = symbols('T', positive=True)
C1a, C2a = Q(-25000), Q(45,10)      # (a) 2H2O(g) <=> 2H2 + O2
C1b, C2b = Q(-31000), Q(155,10)     # (b) 2MO2(s) <=> 2MO(s) + O2
# (c) MO + H2O <=> MO2 + H2 :  (a) = (b) + 2(c)
C1c, C2c = (C1a - C1b)/2, (C2a - C2b)/2
print('(1) log Kc =', C1c, '/T +', C2c)
assert (C1c, C2c) == (3000, Q(-55,10))
R = Q(831,100); k = Q(230,100)*R  # 2.30R [J/(K mol)]
dH = {n: -c*k/1000 for n, c in (('a', C1a), ('b', C1b), ('c', C1c))}
print('(2) ΔH kJ:', {n: N(v, 4) for n, v in dH.items()})
assert abs(dH['b'] + 2*dH['c'] - dH['a']) < 1e-9
# (3) K=1 となる温度
Tb = solve(C1b/T + C2b, T)[0]; Ta = solve(C1a/T + C2a, T)[0]; Tc = solve(C1c/T + C2c, T)[0]
print('(3) Kb=1:', N(Tb), ' Ka=1:', N(Ta), ' Kc=1:', N(Tc, 4))
assert Tb == 2000
# (4) 500 K での Kc と H2 のモル分率（固体過剰，H2O を流す）
lKc = C1c/500 + C2c
Kc = 10**lKc
x = Kc/(1+Kc)
print('(4) 500K log Kc =', lKc, ' Kc =', N(Kc, 4), ' x_H2 =', N(x, 4))
lKc6 = C1c/600 + C2c
print('    600K log Kc =', lKc6, ' x_H2 =', N(10**lKc6/(1+10**lKc6), 3))
assert round(float(x), 2) == 0.76
# (4) 直接熱分解 2000 K での Ka
print('    2000K log Ka =', C1a/2000 + C2a)
# (5) (b) を 1500 K で Kb=1 にする：C2b' = 31000/1500
C2b2 = -C1b/1500
C2c2 = (C2a - C2b2)/2
Tc2 = solve(C1c/T + C2c2, T)[0]
print('(5) C2b\' =', N(C2b2, 4), ' C2c\' =', N(C2c2, 4), ' Kc=1 の温度', N(Tc2, 4))
# 一般式 Tc = 2*3000/(C2b - 4.5)
# 反応式（M に Ce を入れて確認）
for e in ['2H2O(g) <=> 2H2(g) + O2(g)', '2CeO2(s) <=> 2CeO(s) + O2(g)', 'CeO(s) + H2O(g) <=> CeO2(s) + H2(g)']:
    print(e, check_equation(e)); assert check_equation(e) is None
print('OK')
