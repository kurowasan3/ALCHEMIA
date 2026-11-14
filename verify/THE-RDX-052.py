# THE-RDX-052 検証：ポリヨウ化物 Et4N[I(I2)n] の n の決定
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, Rational as R, nsimplify
for e in ['I2 + 2S2O3^2- -> 2I- + S4O6^2-', 'I2 + I- <=> I3-']:
    r = check_equation(e); print(e, r); assert r is None
n = symbols('n', positive=True)
M_salt = 8*12 + 20*1 + 14 + 127          # (C2H5)4NI
assert M_salt == 257
M = M_salt + 254*n
# (2) 0.300 g -> 50.0 mL, 10.0 mL 分取, 0.0200 mol/L で 14.5 mL
nI2 = R('0.0200')*R('14.5')/1000/2
eq = R('0.300')*R(10,50)*n/M - nI2
s = solve(eq, n)[0]; print('滴定 n =', float(s)); assert round(float(s), 1) == 1.6
print('I2 in aliquot', float(nI2))
# (3) 0.500 g -> 残り 0.194 g（Et4NI）
s2 = solve(R('0.500')/M - R('0.194')/M_salt, n)[0]; print('加熱 n =', float(s2)); assert round(float(s2), 1) == 1.6
# (4) I3- (n=1) と I5- (n=2) の混合で平均 1.6 -> I5- の割合 0.6 -> I3-:I5- = 2:3
x = symbols('x'); xs = solve((1-x)*1 + 2*x - R(16,10), x)[0]; print('I5- 割合', xs); assert xs == R(3,5)
# (5) 平均の酸化数 -1/(2n+1)
ox = -1/(2*R(16,10)+1); print('平均酸化数', ox, float(ox)); assert round(float(ox),2) == -0.24
# 理論的に n=1.6 なら滴定量
print('n=1.6 の理論滴定量 mL', float((R('0.0600')*R(16,10)/M.subs(n,R(16,10)))*2/R('0.0200')*1000))
print('OK')
