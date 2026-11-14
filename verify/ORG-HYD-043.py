# 検証スクリプト（東工大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/out/items/ORG-HYD-043/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
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
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
# ---- 個別の検証 ----

import math
r = math.radians
C1 = (0.150*math.cos(r(120)), 0.150*math.sin(r(120)))
C4c = (0.134 + 0.150*math.cos(r(60)), 0.150*math.sin(r(60)))
C4t = (0.134 + 0.150*math.cos(r(-60)), 0.150*math.sin(r(-60)))
d = {'a': math.dist(C1, C4c), 'b': math.dist(C1, C4t), 'c': 0.146*2 + 0.120,
     'd': 0.290, 'e': 0.290*math.sqrt(3), 'f': 0.580}
assert sig(d['a'],3) == 0.284 and sig(d['b'],3) == 0.385 and sig(d['e'],3) == 0.502
order = ''.join(sorted(d, key=lambda k: -d[k])); print(order); assert order == 'fecbda'

print('ALL OK')
