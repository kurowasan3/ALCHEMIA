# POL-PRT-020 検証：ヘキサペプチドの配列と電気泳動
import sys, itertools
sys.path.insert(0, '/home/claude/chem/work/g12')
from g12lib import check_item
CHIRAL = {'Ala': 1, 'Gly': 0, 'Glu': 1, 'Lys': 1, 'Tyr': 1, 'Cys': 1}
def cut(seq, after):
    frags, cur = [], []
    for a in seq:
        cur.append(a)
        if a in after:
            frags.append(cur); cur = []
    if cur: frags.append(cur)
    return frags
sols = []
for p in itertools.permutations(CHIRAL):
    if p[0] != 'Ala' or CHIRAL[p[-1]]:
        continue
    fx = cut(p, {'Lys'})
    if len(fx) != 2: continue
    f2 = fx[1]
    if not (len(f2) <= 2 and 'Cys' in f2): continue
    fy = cut(p, {'Tyr'})
    if sorted(map(len, fy)) != [3, 3]: continue
    sols.append(p)
print(sols)
assert sols == [('Ala', 'Glu', 'Tyr', 'Lys', 'Cys', 'Gly')]
SIDE = {'Glu': -1, 'Lys': +1}
for g in cut(sols[0], {'Tyr'}):
    q = 1 - 1 + sum(SIDE.get(a, 0) for a in g)
    print(g, '正味の電荷', q)
assert [1 - 1 + sum(SIDE.get(a, 0) for a in g) for g in cut(sols[0], {'Tyr'})] == [-1, 1]
assert not check_item('POL-PRT-020')
print('POL-PRT-020 OK')
