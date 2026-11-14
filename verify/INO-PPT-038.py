from fractions import Fraction as F
from itertools import combinations
M1 = {'Na': 23, 'K': 39, 'NH4': 18}
M3 = {'Al': 27, 'Cr': 52, 'Fe': 56}
SO4, H2O, BaSO4 = 96, 18, 137 + 32 + 64
assert BaSO4 == 233
def cand(target):
    out = []
    for a in M1:
        for b in M3:
            if M1[a] + M3[b] == target: out.append(({a: 1}, {b: 1}))
    for a, b in combinations(M1, 2):
        for t in M3:
            y = F(target - M3[t] - M1[b], M1[a] - M1[b])
            if 0 < y < 1: out.append(({a: y, b: 1 - y}, {t: 1}))
    for a, b in combinations(M3, 2):
        for m in M1:
            x = F(target - M1[m] - M3[b], M3[a] - M3[b])
            if 0 < x < 1: out.append(({m: 1}, {a: x, b: 1 - x}))
    return out
data = {'P': F(2265, 10), 'Q': F(2410, 10), 'R': F(2430, 10)}
for k, m in data.items():
    assert (F(1080, 10) / H2O) * 2 == 12        # n = 12
    s = 2*(m - F(108) - SO4)
    data[k] = s
print({k: int(v) for k, v in data.items()})
def flt(c, nh4, flame_ok, fe):
    return [x for x in c if (('NH4' in x[0]) == nh4) and (set(x[0]) - {'NH4'} <= flame_ok) and (('Fe' in x[1]) == fe)]
P = flt(cand(data['P']), True, set(), False)
Q = flt(cand(data['Q']), True, set(), True)
print('P', P); print('Q', Q)
assert P == [({'NH4': 1}, {'Al': 1})] and Q == [({'NH4': 1}, {'Fe': 1})]
# R: (i) NH4 なし (ii) K のみ (iii) Fe なし
R_ii = [x for x in cand(data['R']) if set(x[0]) == {'K'}]
print('R, K only', R_ii)
assert len(R_ii) == 2
R = [x for x in R_ii if 'Fe' not in x[1]]
assert R == [({'K': 1}, {'Al': F(13, 25), 'Cr': F(12, 25)})]
R4 = [x for x in cand(data['R']) if 'NH4' not in x[0] and 'Fe' not in x[1] and x not in R]
print('without (ii)', R4)
assert R4 == [({'Na': F(13, 16), 'K': F(3, 16)}, {'Cr': 1})]
Al2O3, Cr2O3 = 2*27 + 48, 2*52 + 48
m3 = F(1, 2)    # mmol
ox3 = m3*F(13, 25)/2*Al2O3 + m3*F(12, 25)/2*Cr2O3
ox4 = m3/2*Cr2O3
print('oxide mg', float(ox3), float(ox4))
assert round(float(ox3), 1) == 31.5 and float(ox4) == 38.0
print('OK')
