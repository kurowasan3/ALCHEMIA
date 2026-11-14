# THE-BND-097 検証：KBr ナノ結晶（3x3 断面）のイオン数・質量・包む円筒の最小直径
from sympy import Rational, sqrt, N
def count(nx, ny, nz):
    K = sum(1 for i in range(nx) for j in range(ny) for k in range(nz) if (i+j+k) % 2 == 0)
    return K, nx*ny*nz - K
for n in range(1, 9):
    K, B = count(3, 3, n)
    print('3x3x%d' % n, K, B, 'neutral' if K == B else 'charge %+d' % (K-B))
for n in range(1, 5):
    print('4x4x%d' % n, count(4, 4, n))
# (2) 3x3x6 の質量
K, B = count(3, 3, 6); assert K == B == 27
M = 39.1 + 79.9
m = 27*M/6.02e23
print('mass g', m)               # 5.337e-21 -> 5.34e-21
assert f'{m:.2e}' == '5.34e-21'
# (3) 最小直径：中心のイオンを円の中心にとる。角（対角線上）は中心から d*sqrt2，辺の中央は d
rK, rBr = 0.138, 0.196
d = rK + rBr
s2 = 1.414
corner_Br = d*s2 + rBr   # 角が Br- の層
corner_K = d*s2 + rK
edge_K = d + rK          # 角が Br- の層では辺の中央は K+
edge_Br = d + rBr
R = max(corner_Br, corner_K, edge_K, edge_Br)
print('R', R, 'D', 2*R)
assert f'{2*R:.3g}' == '1.34'
# sqrt2 を正確に使っても 1.34
print('exact', N(2*(d*sqrt(2)+rBr)))
