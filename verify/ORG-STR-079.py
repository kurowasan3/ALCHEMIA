# ORG-STR-079 verification: ESR hyperfine patterns by convolution
from math import comb
from collections import defaultdict
def spectrum(sets):
    """sets: list of (n equivalent H, coupling in units of 0.1 mT)"""
    lines = {0: 1}
    for n, a in sets:
        new = defaultdict(int)
        for p, I in lines.items():
            for k in range(n + 1):
                new[p + (k - n / 2) * a] += I * comb(n, k)
        lines = dict(new)
    return sorted(lines.items())
aA, aB = 22, 25
cases = {
    'methyl': [(3, 23)],
    'isopropyl': [(1, aA), (6, aB)],
    'n-propyl': [(2, aA), (2, aB)],
    'tert-butyl': [(9, aB)],
    'isobutyl': [(2, aA), (1, aB)],
    'cyclopentyl': [(1, aA), (4, aB)],
    'isopropyl_equal': [(1, aB), (6, aB)],
}
R = {}
for k, s in cases.items():
    sp = spectrum(s); I = [i for p, i in sp]
    R[k] = (len(sp), max(I) // min(I), I)
    print(k, len(sp), max(I) / min(I), I)
assert R['methyl'][:2] == (4, 3)
assert R['isopropyl'][:2] == (14, 20) and R['n-propyl'][0] == 9
assert R['tert-butyl'][:2] == (10, 126) and R['isobutyl'][0] == 6
assert R['cyclopentyl'][:2] == (10, 6)
assert R['isopropyl_equal'][0] == 8 and R['isopropyl_equal'][2] == [1, 7, 21, 35, 35, 21, 7, 1]
# equal spacing check for tert-butyl vs cyclopentyl
pt = [p for p, i in spectrum(cases['tert-butyl'])]; pc = [p for p, i in spectrum(cases['cyclopentyl'])]
dt = {round(b - a, 6) for a, b in zip(pt, pt[1:])}; dc = sorted({round(b - a, 6) for a, b in zip(pc, pc[1:])})
print('spacings tBu', dt, 'cyclopentyl', dc)
assert dt == {25.0} and dc == [3.0, 22.0]
# molecular formulas: isobutane C4H10, cyclopentane C5H10; radicals C4H9, C5H9
print('OK')
