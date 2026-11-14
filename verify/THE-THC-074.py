# THE-THC-074 verification: cyclic sulfur imides, isomer enumeration and bond-energy enthalpy
from itertools import combinations
from collections import Counter
E = dict(SS=226, SN=240, NN=160, NH=391, HH=436, N2=946)
def canon(s):
    v = []
    for r in range(8):
        t = s[r:] + s[:r]; v += [t, t[::-1]]
    return min(v)
def dH(s):
    n = sum(s)
    bonds = Counter()
    for i in range(8):
        a, b = s[i], s[(i + 1) % 8]
        bonds['NN' if a and b else ('SN' if a or b else 'SS')] += 1
    broken = (8 - n) * E['SS'] + n / 2 * E['N2'] + n / 2 * E['HH']
    formed = sum(bonds[k] * E[k] for k in bonds) + n * E['NH']
    return broken - formed, bonds['NN']
res = {}
for n in (2, 3, 4):
    iso = {}
    for c in combinations(range(8), n):
        s = tuple(1 if i in c else 0 for i in range(8))
        iso[canon(s)] = dH(s)
    for s, (h, k) in iso.items():
        assert abs(h - (46 * n + 94 * k)) < 1e-9
    res[n] = (len(iso), Counter(k for h, k in iso.values()), sorted(h for h, k in iso.values()))
    print(n, res[n])
assert res[2][0] == 4 and res[2][1] == Counter({0: 3, 1: 1})
assert res[3][0] == 5 and res[3][1] == Counter({0: 2, 1: 2, 2: 1})
assert res[4][0] == 8 and res[4][1] == Counter({0: 1, 1: 2, 2: 4, 3: 1})
assert res[4][2][0] == 184 and res[4][2][-1] == 466
print('OK')
