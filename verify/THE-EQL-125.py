# THE-EQL-125 検証：ジ-tert-ブチルペルオキシドの一次分解と全圧
import math, sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
assert check_equation('C8H18O2 -> 2C3H6O + C2H6') is None
# P = P0(1 + 2x)
x20 = (1.50 - 1)/2; r20 = 1 - x20
print('remaining at 20 min', r20)
r40 = r20**2; P40 = 1 + 2*(1 - r40)
print('remaining 40', r40, 'P40/P0', P40)
assert abs(P40 - 1.875) < 1e-12
log2, log3 = 0.301, 0.477
logr = log3 - 2*log2
t_half = 20*log2/(-logr)
print('t1/2 (given logs)', t_half, 'exact', 20*math.log(2)/(-math.log(0.75)))
assert round(t_half) == 48
# 20-40 分の平均速度の比（アセトン生成/過酸化物分解）= 2
print('ratio 2 (stoichiometry)')
# 初圧を 2 倍にしても全圧が 1.50 倍になる時間は 20 分（一次反応）
print('OK')
