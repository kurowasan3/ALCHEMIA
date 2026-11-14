# 検証スクリプト THE-BND-033  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-BND-033/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b:
                eqs.add(b)
    return sorted(eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
import itertools
import numpy as np
from sympy import sqrt, pi, Rational as R, N
check_structure()
assert all_equations() == []
r = 1.0
a = 4 / np.sqrt(3)
atoms = [np.array(p, float) * a for p in itertools.product(range(-2, 4), repeat=3)]
atoms += [(np.array(p, float) + 0.5) * a for p in itertools.product(range(-2, 3), repeat=3)]
def shells(site):
    s = np.array(site) * a
    d = sorted(round(float(np.linalg.norm(x - s)), 4) for x in atoms)
    out = {}
    for x in d[:12]: out[x] = out.get(x, 0) + 1
    return out
O = shells((0.5, 1.0, 0.5)); T = shells((0.5, 1.0, 0.25)); Oe = shells((0.5, 1.0, 0.0))
print('O', O); print('T', T); print('O(辺)', Oe)
assert list(O.items())[:2] == [(1.1547, 2), (1.633, 4)] and O == Oe or list(Oe.items())[:2] == [(1.1547, 2), (1.633, 4)]
assert list(T.items())[0] == (1.291, 4) and list(T.items())[1][0] == 2.0817
# すき間の個数（単位格子内の位置を周期境界で数える）
def count_sites(cands):
    pts = set()
    for c in cands:
        for perm in set(itertools.permutations(c)):
            pts.add(tuple(x % 1 for x in perm))
    return pts
Osites = count_sites([(0.5, 0.5, 0.0), (0.5, 0.0, 0.0)])
Tsites = count_sites([(0.5, 0.25, 0.0), (0.5, 0.75, 0.0)])
print(len(Osites), len(Tsites)); assert len(Osites) == 6 and len(Tsites) == 12
# 半径
ab = 4 / sqrt(3)
oO = ab/2 - 1; oT = sqrt(5)*ab/4 - 1; af = 2*sqrt(2); fo = af/2 - 1; ft = sqrt(3)*af/4 - 1
assert sig(float(oO), 3) == 0.155 and sig(float(oT), 3) == 0.291 and sig(float(fo), 3) == 0.414 and sig(float(ft), 3) == 0.225
assert sig(2/1.732 - 1, 3) == 0.155 and sig(3.873/3 - 1, 3) == 0.291
phib = 1.732*3.14/8; phif = 0.740; print('φb', phib); assert sig(phib, 3) == 0.680
assert sig(float(sqrt(2)*pi/6), 3) == 0.740
vb = (1/0.680 - 1) / 9; vf = (1/0.740 - 1) / 3; print(vb, vf)
assert sig(vb, 2) == 0.052 and sig(vf, 2) == 0.12
rFe, rC = 0.124, 0.077
ratios = [rC/(0.414*rFe), rC/(0.291*rFe), rC/(0.155*rFe)]; print(ratios)
assert [sig(x, 3) for x in ratios] == [1.50, 2.13, 4.01]
assert [sig(rC/(float(h)*rFe), 3) for h in (fo, oT, oO)] == [1.50, 2.13, 4.01]
occ = (2.1/12) / (97.9/56) * 100; print('occ', occ); assert sig(occ, 2) == 10
print('全検証OK:', ITEM)
