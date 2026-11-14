# THE-BND-110 verification: bcc iron surfaces, coordination numbers, model indices
import itertools, numpy as np
from math import sqrt
R = 6
P = []
for i, j, k in itertools.product(range(-R, R + 1), repeat=3):
    P.append((i, j, k)); P.append((i + .5, j + .5, k + .5))
P = np.array(P)
NN = [np.array(v) * 0.5 for v in itertools.product((1, -1), repeat=3)]
res = {}
for h in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 1)]:
    n = np.array(h, float); n /= np.linalg.norm(n)
    d = 1 / np.linalg.norm(h) if h != (1, 0, 0) else 0.5   # interlayer spacing (a=1) for bcc
    if h == (1, 1, 0): d = 1 / sqrt(2)
    if h == (1, 1, 1): d = 1 / (2 * sqrt(3))
    if h == (2, 1, 1): d = 1 / sqrt(6)
    # layer distribution of NN vectors
    dist = sorted(round(float(v @ n) / d) for v in NN)
    depth = -(P @ n); keep = depth > -1e-9
    S = set(map(tuple, np.round(P[keep] * 2).astype(int)))
    layers = {}
    for p, dd in zip(P[keep], depth[keep]):
        if dd > 1.6 or np.max(np.abs(p)) > R - 2: continue
        cn = sum(1 for v in NN if tuple(np.round((p + v) * 2).astype(int)) in S)
        if cn < 8: layers.setdefault(int(round(dd / d)) + 1, set()).add(cn)
    sigma = 2 * d
    missing = sum(8 - list(c)[0] for c in layers.values())
    c7 = sum(1 for c in layers.values() if 7 in c)
    res[h] = (dist, layers, round(sigma, 3), round(missing * sigma, 2), round(c7 * sigma, 3))
    print(h, res[h])
assert res[(1, 0, 0)][1] == {1: {4}} and res[(1, 1, 0)][1] == {1: {6}}
assert res[(1, 1, 1)][1] == {1: {4}, 2: {7}, 3: {7}} and res[(2, 1, 1)][1] == {1: {5}, 2: {7}}
assert [res[h][3] for h in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 1)]] == [4.0, 2.83, 3.46, 3.27]
assert [res[h][4] for h in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 1)]] == [0, 0, 1.155, 0.816]
assert res[(1, 1, 1)][0] == [-3, -1, -1, -1, 1, 1, 1, 3]
assert res[(2, 1, 1)][0] == [-2, -1, -1, 0, 0, 1, 1, 2]
assert res[(1, 1, 0)][0] == [-1, -1, 0, 0, 0, 0, 1, 1]
a = 0.287e-7
print('110 per cm2', sqrt(2) / a**2); assert round(sqrt(2) / a**2 / 1e15, 1) == 1.7
print('OK')
