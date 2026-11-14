# 検証スクリプト（東工大 g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/out/items/THE-BND-079/item.json'
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans <= 1), (nprob, nans)
print('小問数', nprob, '答え', nans)
for k in ['src_no', 'src_label', 'field', 'theme', 'difficulty', 'time_min', 'format', 'problem', 'answer', 'solution', 'core', 'tags', 'verified', 'verify_note', 'problem_fig', 'solution_fig']:
    assert k in it, k
# ---- 個別の検証 ----

import itertools, math
R = range(-4, 5)
# NaCl: half-a units; site (i,j,k) in units a/2: Cl if i+j+k even, Na if odd
def d(p): return math.sqrt(sum(c*c for c in p))
cl = [p for p in itertools.product(R, R, R) if sum(p) % 2 == 0 and p != (0, 0, 0)]
na = [p for p in itertools.product(R, R, R) if sum(p) % 2 == 1]
from collections import Counter
cc = Counter(round(d(p)/2, 4) for p in cl); cn = Counter(round(d(p)/2, 4) for p in na)
print(sorted(cc.items())[:2], sorted(cn.items())[:2])
assert sorted(cc.items())[0] == (round(math.sqrt(2)/2, 4), 12)
assert sorted(cn.items())[0] == (0.5, 6) and sorted(cn.items())[1] == (round(math.sqrt(3)/2, 4), 8)
# CsCl: Cl at integer points, Cs at half-integer (unit b)
clc = [p for p in itertools.product(R, R, R) if p != (0, 0, 0)]
cs = [tuple(c + 0.5 for c in p) for p in itertools.product(R, R, R)]
c1 = Counter(round(d(p), 4) for p in clc); c2 = Counter(round(d(p), 4) for p in cs)
assert sorted(c1.items())[0] == (1.0, 6) and sorted(c2.items())[0] == (round(math.sqrt(3)/2, 4), 8)
near_cs = [p for p in cs if abs(d(p) - math.sqrt(3)/2) < 1e-6]
S = set()
for c in near_cs:
    for q in itertools.product(R, R, R):
        if abs(d(tuple(q[i] - c[i] for i in range(3))) - math.sqrt(3)/2) < 1e-6 and q != (0, 0, 0): S.add(q)
print('CsCl count', len(S)); assert len(S) == 26
near_na = [p for p in na if abs(d(p)/2 - 0.5) < 1e-6]
S2 = set()
for c in near_na:
    for q in itertools.product(R, R, R):
        if sum(q) % 2 == 0 and q != (0, 0, 0) and abs(d(tuple(q[i] - c[i] for i in range(3)))/2 - 0.5) < 1e-6: S2.add(q)
print('NaCl count', len(S2)); assert len(S2) == 18
assert round(1.41/2, 3) == 0.705 and round(1.73/2, 3) == 0.865
print(round(math.sqrt(2), 2), round(2/math.sqrt(3), 2)); assert round(2/math.sqrt(3), 2) == 1.15
print('OK')
