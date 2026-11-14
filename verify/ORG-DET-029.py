"""ORG-DET-029 検証：陽性表・方式1で必要な試験・方式2の最悪回数の最小値・示した手順"""
import sys, json, itertools, re
from functools import lru_cache
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
T = ['I', 'II', 'III', 'IV', 'V']   # Na, NaHCO3, 銀鏡, ヨードホルム, FeCl3
D = {'エタノール': '10010', 'アセトアルデヒド': '00110', 'アセトン': '00010', '酢酸': '11000',
     'ギ酸': '11100', 'フェノール': '10001', 'メタノール': '10000', 'ジエチルエーテル': '00000'}
# tests.json と照合
tj = {t['id']: t for t in json.load(open('/home/claude/kurowasan3/alchemia/content/chem/tests.json'))['tests']}
idmap = {'エタノール': 'ethanol', 'アセトアルデヒド': 'acetaldehyde', 'アセトン': 'acetone', '酢酸': 'acetic-acid',
         'ギ酸': 'formic-acid', 'フェノール': 'phenol', 'メタノール': 'methanol', 'ジエチルエーテル': 'diethyl-ether'}
tid = ['sodium', 'nahco3', 'silver-mirror', 'iodoform', 'fecl3']
for name, v in D.items():
    for i, t in enumerate(tid):
        cid = idmap[name]
        if cid in tj[t]['scope']:
            assert (cid in tj[t]['pos']) == (v[i] == '1'), (name, t)
        else:
            print('  scope外:', name, t, v[i])
names = list(D)
# 方式1
sep = [s for k in range(1, 6) for s in itertools.combinations(range(5), k)
       if len({tuple(D[n][i] for i in s) for n in names}) == len(names)]
print('区別できる組:', sep); assert sep == [(0, 1, 2, 3, 4)]
for a, b in itertools.combinations(names, 2):
    diff = [T[i] for i in range(5) if D[a][i] != D[b][i]]
    if len(diff) == 1: print('  1か所だけ異なる:', a, b, diff)
# 方式2
@lru_cache(None)
def depth(S, avail):
    if len(S) <= 1: return 0
    best = 99
    for t in avail:
        g1 = frozenset(c for c in S if D[c][t] == '1'); g0 = S - g1
        if g0 and g1:
            rest = tuple(x for x in avail if x != t)
            best = min(best, 1 + max(depth(g0, rest), depth(g1, rest)))
    return best
d = depth(frozenset(names), tuple(range(5))); print('方式2 最小の最悪回数', d); assert d == 4
print('陽性数', [sum(D[n][i] == '1' for n in names) for i in range(5)])
# 示した手順
def proc(v):
    seq = [0]
    if v[0] == '0':
        seq.append(3)
        if v[3] == '1': seq.append(2)
    else:
        seq.append(1)
        if v[1] == '1': seq.append(2)
        else:
            seq.append(3)
            if v[3] == '0': seq.append(4)
    return seq
res = {}
for n in names:
    s = proc(D[n]); key = tuple((i, D[n][i]) for i in s); res.setdefault(key, []).append(n)
    assert len(s) <= 4
assert all(len(v) == 1 for v in res.values()); print('手順で全8種類が決まる')
assert '10010' == D['エタノール']  # 2-プロパノールも Na+, ヨードホルム+ → 10010
for eq in ['CH3CH2OH -> CH3CHO + H2', '2CH3CH(OH)CH3 + O2 -> 2CH3COCH3 + 2H2O']:
    assert check_equation(eq) is None
print('ALL OK')
