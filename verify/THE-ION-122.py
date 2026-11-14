# THE-ION-122 検証：固体ヨウ素と KI 水溶液（I3- 生成）の溶解平衡
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, nsolve, Rational as R, N as Nn
assert check_equation('I2 + I- <=> I3-') is None
K1 = R('1.3e-3')
# (1) 0.100 mol/L KI の飽和溶液の全ヨウ素（I2 換算）4.89e-2 mol/L -> K2
S = R('4.89e-2'); c = R('0.100')
I3 = S - K1; Im = c - I3
K2 = I3/(K1*Im); print('K2 =', float(K2)); assert round(float(K2), -1) == 700
K2 = R(700)
# (2) 0.500 mol/L KI の溶解度と S/c の上限
def sat(c):
    Im = c/(1+K1*K2); I3 = K1*K2*Im; return K1 + I3, Im, I3
S05 = sat(R('0.500'))[0]; print('S(0.500) =', float(S05)); assert round(float(S05), 3) == 0.240
lim = K1*K2/(1+K1*K2); print('上限比', float(lim)); assert round(float(lim), 2) == 0.48
# (3) 0.0200 mol/L KI 1.00 L: 固体が消える境目 a*
c = R('0.0200')
Sst, Im_s, I3_s = sat(c); print('a* =', float(Sst), '[I-]=', float(Im_s), '[I3-]=', float(I3_s))
assert round(float(Sst), 4) == 0.0108
# (3a) a = 5.00e-3 mol -> 全部溶ける：x=[I3-], K2 (a-x)(c-x) = x
a = R('5.00e-3'); x = symbols('x')
sols = [s for s in solve(K2*(a-x)*(c-x) - x, x) if 0 < s < a]
xs = sols[0]; I2aq = a - xs; Imin = c - xs
print('x=[I3-]', float(xs), '[I2]', float(I2aq), '[I-]', float(Imin))
assert I2aq < K1   # 固体が残らない条件の確認
# 近似（x≈a）の [I2]
I2approx = a/(K2*(c-a)); print('近似 [I2]', float(I2approx), '相対誤差', float((I2approx-I2aq)/I2aq))
# (3b) a = 2.00e-2 mol -> 固体が残る： 残る固体の質量
rest = (R('2.00e-2') - Sst)*254; print('残る固体 g', float(rest)); assert round(float(rest), 2) == 2.33
# (4) KI 1 mol あたり溶ける I2 の量の上限：c が大きいとき S/c -> 0.476
for cc in [R('0.0200'), R('0.100'), R('0.500'), R('1.00')]:
    s_, *_ = sat(cc); print(float(cc), float(s_/cc))
print('OK')
