"""THE-BND-067 検証"""
import itertools
from fractions import Fraction as F
from math import gcd
# 電荷のつり合い
def ratio(zm, zc):
    a, b = zc, zm   # zm*a = zc*b
    g = gcd(a, b); return a // g, b // g
assert ratio(3, 4) == (4, 3) and ratio(2, 3) == (3, 2)
assert F(4 - 3, 4) == F(1, 4) and F(3 - 2, 3) == F(1, 3)
# NaCl 型：C 側（i+j+k 奇数）の最近接 N 側の数，単位格子中の各位置の数
pts = list(itertools.product(range(-1, 4), repeat=3))
c = (1, 0, 0)
nn = [p for p in pts if sum((p[i] - c[i]) ** 2 for i in range(3)) == 1]
assert len(nn) == 6 and all(sum(p) % 2 == 0 for p in nn)
cnt = {0: F(0), 1: F(0)}
for i, j, k in itertools.product(range(3), repeat=3):
    w = F(1)
    for t in (i, j, k):
        if t in (0, 2): w /= 2
    cnt[(i + j + k) % 2] += w
assert cnt == {0: 4, 1: 4}
M1 = 55.8 * 4 + (55.8 + 26.0 * 6) * 3 + 18.0 * 6
M2 = 39.1 + 55.8 * 2 + 26.0 * 6
M3 = 63.5 * 3 + (55.8 + 156.0) * 2 + 108.0
V = (1.02e-7) ** 3
r1, r2 = M1 / 6.02e23 / V, 4 * M2 / 6.02e23 / V
loss = 108.0 / M3 * 100
print(M1, M2, M3, r1, r2, loss)
assert abs(M1 - 966.6) < 1e-6 and abs(M2 - 306.7) < 1e-6 and abs(M3 - 722.1) < 1e-6
assert f'{r1:.2g}' == '1.5' and f'{r2:.2g}' == '1.9' and f'{loss:.2g}' == '15'
assert 18 / 4 == 4.5
print('THE-BND-067 OK')
