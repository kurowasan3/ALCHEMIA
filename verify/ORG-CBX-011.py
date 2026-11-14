"""ORG-CBX-011 検証：油脂の分子量・C=C 数・構成脂肪酸・不斉炭素"""
import sys, itertools
sys.path.insert(0, '/home/claude/chem/work/g09')
from sympy import Rational, nsimplify
from g09lib import report_equations, formula, n_stereocenters
from rdkit.Chem.Descriptors import MolWt

M = {'H': Rational(1), 'C': 12, 'O': 16, 'K': 39}
nKOH = Rational(168, 100) / (M['K'] + M['O'] + M['H'])
nX = nKOH / 3
MX = Rational(856, 100) / nX
print('nKOH', nKOH, 'MX', MX); assert MX == 856
nH2 = Rational(672, 1000) / Rational(224, 10)
print('C=C per molecule', nH2 / nX); assert nH2 / nX == 3

acids = {'palmitic': (16, 0), 'stearic': (18, 0), 'oleic': (18, 1), 'linoleic': (18, 2), 'linolenic': (18, 3)}
def mw(c, db):  # CnH(2n-2db)O2
    return 12 * c + (2 * c - 2 * db) + 32
for k, (c, db) in acids.items():
    print(k, mw(c, db))
glyc = 92
target = MX - glyc + 54 - mw(16, 0)
pairs = [(a, b) for a, b in itertools.combinations_with_replacement([k for k in acids if acids[k][0] == 18], 2)
         if mw(*acids[a]) + mw(*acids[b]) == target and acids[a][1] + acids[b][1] == 3]
print('C18 pairs summing to', target, 'with 3 C=C:', pairs)
assert set(pairs) == {('stearic', 'linolenic'), ('oleic', 'linoleic')}
ok = [p for p in pairs if 'stearic' not in p]   # 飽和脂肪酸はパルミチン酸だけ
assert ok == [('oleic', 'linoleic')]

acyl = {'P': 'C(=O)CCCCCCCCCCCCCCC', 'S': 'C(=O)CCCCCCCCCCCCCCCCC',
        'O': 'C(=O)CCCCCCCC=CCCCCCCCC', 'L': 'C(=O)CCCCCCCC=CCC=CCCCCC'}
for k, a in acyl.items():
    print(k, formula('O' + a))
def tg(e1, mid, e2):
    return f'C(O{acyl[e1]})C(O{acyl[mid]})CO{acyl[e2]}'
xs = {'P-mid': tg('O', 'P', 'L'), 'O-mid': tg('P', 'O', 'L'), 'L-mid': tg('P', 'L', 'O')}
hyd = {'P-mid': tg('S', 'P', 'S'), 'O-mid': tg('P', 'S', 'S'), 'L-mid': tg('P', 'S', 'S')}
from rdkit import Chem
for k in xs:
    m = Chem.MolFromSmiles(xs[k])
    print(k, 'X', formula(xs[k]), round(MolWt(m), 1), 'stereo', n_stereocenters(xs[k]),
          '| Y', formula(hyd[k]), 'stereo', n_stereocenters(hyd[k]))
    assert formula(xs[k]) == 'C55H100O6' and n_stereocenters(xs[k]) == 1
assert n_stereocenters(hyd['P-mid']) == 0 and n_stereocenters(hyd['O-mid']) == 1
# 整数原子量での分子量
print('C55H100O6 =', 55 * 12 + 100 + 96)
assert 55 * 12 + 100 + 96 == 856
report_equations('ORG-CBX-011')
