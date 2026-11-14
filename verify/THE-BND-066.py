# THE-BND-066 検証：NaCl 型結晶のマーデルング定数（殻ごとの和と電気的に中性な立方体による和）
import itertools, math
from sympy import sqrt, Rational as Q, nsimplify, N
# 中心の陽イオンから見た位置 (i,j,k)（単位は最近接距離 r）、電荷は (-1)^(i+j+k)
def shells(nmax):
    d = {}
    for i, j, k in itertools.product(range(-nmax, nmax+1), repeat=3):
        if (i, j, k) == (0, 0, 0): continue
        n2 = i*i + j*j + k*k
        d.setdefault(n2, [0, (-1)**((i+j+k) % 2)])[0] += 1
    return d
sh = shells(4)
first = sorted(sh)[:4]
print('shells (dist^2, number, sign):', [(n2, sh[n2][0], sh[n2][1]) for n2 in first])
assert [(n2, sh[n2][0], sh[n2][1]) for n2 in first] == [(1, 6, -1), (2, 12, 1), (3, 8, -1), (4, 6, 1)]
# 部分和 A_m = -Σ sign*N/sqrt(n2)（引力側を正）
A = 0; partial = []; charge = 1; charges = []
for n2 in first:
    N_, s = sh[n2]; A += -s*N_/sqrt(n2); partial.append(A); charge += s*N_; charges.append(charge)
vals = [round(float(a), 3) for a in partial]
print('partial sums:', vals, ' net charge inside:', charges)
assert vals == [6.0, -2.485, 2.134, -0.866] and charges == [-5, 7, -1, 5]
# エフィェン（Evjen）の立方体：辺 2r（27 個）と 辺 4r（125 個）
def evjen(m):
    A = 0; q = 0
    for i, j, k in itertools.product(range(-m, m+1), repeat=3):
        w = Q(1)
        for c in (i, j, k):
            if abs(c) == m: w *= Q(1, 2)
        s = (-1)**((i+j+k) % 2)
        q += s*w
        if (i, j, k) != (0, 0, 0):
            A += -s*w/sqrt(i*i + j*j + k*k)
    return A, q
A1, q1 = evjen(1); A2, q2 = evjen(2)
print('cube 2r: A =', nsimplify(A1), '=', N(A1, 5), ' charge', q1)
print('cube 4r: A =', N(A2, 6), ' charge', q2)
assert q1 == 0 and q2 == 0 and round(float(A1), 3) == 1.456 and round(float(A2), 2) == 1.75
# 辺 4r の立方体の殻ごとの内訳（解説用）
m = 2; tab = {}
for i, j, k in itertools.product(range(-m, m+1), repeat=3):
    if (i, j, k) == (0, 0, 0): continue
    key = tuple(sorted((abs(i), abs(j), abs(k)), reverse=True))
    w = Q(1)
    for c in (i, j, k):
        if abs(c) == m: w *= Q(1, 2)
    tab.setdefault(key, [0, w, (-1)**((i+j+k) % 2)])[0] += 1
tot = 0
for key in sorted(tab, key=lambda t: sum(c*c for c in t)):
    n, w, s = tab[key]; term = -s*n*w/sqrt(sum(c*c for c in key)); tot += term
    print(key, 'N=', n, 'w=', w, 'sign', s, 'term', round(float(term), 4))
print('sum =', round(float(tot), 4)); assert abs(float(tot) - float(A2)) < 1e-12
# 文献値 1.7476 との比較
print('known Madelung constant NaCl = 1.7476')
print('OK')
