# POL-PRT-019 検証：ヘキサペプチドの配列
import sys, itertools
sys.path.insert(0, '/home/claude/chem/work/g12')
from g12lib import check_item
def trypsin(seq):
    frags, cur = [], []
    for a in seq:
        cur.append(a)
        if a == 'Lys':
            frags.append(cur); cur = []
    if cur: frags.append(cur)
    return frags
sols = []
for p in itertools.permutations(['Gly', 'Ala', 'Ser', 'Glu', 'Lys', 'Phe']):
    if p[0] != 'Glu' or p[-3:] != ('Ala', 'Ser', 'Phe'):
        continue
    if sorted(len(f) for f in trypsin(p)) == [3, 3]:
        sols.append(p)
print(sols)
assert sols == [('Glu', 'Gly', 'Lys', 'Ala', 'Ser', 'Phe')]
alt = trypsin(['Glu', 'Lys', 'Gly', 'Ala', 'Ser', 'Phe']); print(alt)
assert alt == [['Glu', 'Lys'], ['Gly', 'Ala', 'Ser', 'Phe']]
biuret = [f for f in alt if len(f) - 1 >= 2]; assert biuret == [['Gly', 'Ala', 'Ser', 'Phe']]
assert not check_item('POL-PRT-019')
print('POL-PRT-019 OK')
