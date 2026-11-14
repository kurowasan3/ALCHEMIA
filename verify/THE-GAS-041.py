# THE-GAS-041 検証：理想気体のグラフの読み取り（n = 0.20 mol, R = 8.3e3）
from sympy import Rational as Q, N
n, R = Q(20,100), Q(83,10)*1000
# (a) 400 K で V = 6, 4, 2 L の直線が p2, p3, p1
V = {'p2':6, 'p3':4, 'p1':2}
p = {k: n*R*400/v for k, v in V.items()}
print('(a)', {k: N(v,4) for k, v in p.items()})
order = sorted(p, key=p.get); print('   小さい順', order); assert order == ['p2','p3','p1']
r = [p['p1'], p['p2'], p['p3']]; g = min(r); print('   比', [x/g for x in r])
assert [x/g for x in r] == [3, 1, Q(3,2)]  # → 6:2:3
assert abs(float(p['p1']) - 3.3e5) < 0.05e5
# (b) 1/V = 0.50 L^-1 で p = 4, 2, 3 ×10^5 Pa が TA, TB, TC
P = {'TA':4e5, 'TB':2e5, 'TC':3e5}
T = {k: v*2.0/float(n*R) for k, v in P.items()}  # V = 2.0 L
print('(b)', {k: round(v,1) for k, v in T.items()})
assert sorted(T, key=T.get) == ['TB','TC','TA'] and round(T['TA'], -1) == 480
# (3) pV/T は温度によらず nR（同じ水平線）
print('(3) pV/T = nR =', N(n*R,4))
# (4) V = (nR/p)(t + 273.15) → t 切片 -273 ℃
print('(4) t 切片 = -273 ℃')
print('OK')
