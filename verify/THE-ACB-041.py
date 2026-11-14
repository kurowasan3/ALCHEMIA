# THE-ACB-041 ナノ空間の pH：独立計算
from sympy import Rational as R, pi, log, N, sqrt
NA = 6.02e23
pi_ = 3.14
d = 6.0e-9            # m 内径
r = d/2
V_m3 = 4/3*pi_*r**3
V_L = V_m3*1e3
print('V (L)=', V_L)                       # 1.13e-22
c1 = 1/(NA*V_L); print('c(1 H3O+) =', c1, 'pH =', -float(N(log(c1,10))))
# (2) pH 5.0 の外液と平衡：カプセル1個あたりの平均 H3O+ 数
avg = 1e-5*NA*V_L; print('avg per capsule =', avg, '-> 1 in', 1/avg)
# (3) 内径を半分
V2 = V_L/8; c2 = 1/(NA*V2); print('half d: pH =', -float(N(log(c2,10))), 'ΔpH =', float(N(log(8,10))))
avg2 = 1e-5*NA*V2; print('half d: 1 in', 1/avg2)
# 純水 pH7 では
print('pure water avg', 1e-7*NA*V_L, '1 in', 1/(1e-7*NA*V_L))
