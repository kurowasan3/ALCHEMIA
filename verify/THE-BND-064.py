"""THE-BND-064 検証：格子点の数え上げ・電荷・密度"""
from fractions import Fraction as F
import itertools
# 一辺 1.02 nm（=2 単位）の立方体，格子点 (i,j,k)∈{0,1,2}^3，(i+j+k) 偶数=Ni，奇数=Fe
cnt = {'Ni': F(0), 'Fe': F(0)}
for i, j, k in itertools.product(range(3), repeat=3):
    w = F(1)
    for t in (i, j, k):
        if t in (0, 2): w /= 2
    cnt['Ni' if (i + j + k) % 2 == 0 else 'Fe'] += w
assert cnt == {'Ni': 4, 'Fe': 4}, cnt
for fe, x in ((2, 2), (3, 1)):
    assert fe - 6 + 2 + x == 0
M = 39.1 * 2 + 58.7 + 55.8 + (12.0 + 14.0) * 6
rho = 4 * M / 6.02e23 / (1.02e-7) ** 3
print(M, rho)
assert abs(M - 348.7) < 1e-9 and f'{rho:.2g}' == '2.2'
assert 2 * 4 / 8 == 1.0 and 1 * 4 / 8 == 0.5
print('THE-BND-064 OK')
