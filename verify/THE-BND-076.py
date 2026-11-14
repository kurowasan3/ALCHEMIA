# 検証スクリプト（東工大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/out/items/THE-BND-076/item.json'
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

import itertools, math
B = [p for p in itertools.product((0,1),repeat=3)] + [(.5,.5,.5)]
A = [(.25,.25,.25),(.75,.75,.25),(.75,.25,.75),(.25,.75,.75)]
def images(p):
    for d in itertools.product((-1,0,1),repeat=3):
        yield tuple(pi+di for pi,di in zip(p,d))
allB = set(q for p in B for q in images(p))
dA = [math.dist(A[0], q) for q in allB]; m = min(dA)
assert abs(m - math.sqrt(3)/4) < 1e-9 and sum(1 for d in dA if abs(d-m) < 1e-9) == 2
allA = set(q for p in A for q in images(p))
dB = [math.dist((.5,.5,.5), q) for q in allA]; assert sum(1 for d in dB if abs(d-m) < 1e-9) == 4
dB0 = [math.dist((0,0,0), q) for q in allA]; assert sum(1 for d in dB0 if abs(d-m) < 1e-9) == 4
assert abs(min(math.dist(A[0],q) for q in allA if q != A[0]) - 1/math.sqrt(2)) < 1e-9
assert round(1.73*0.472/4,3) == 0.204
assert sig(2*232/(6.0e23*0.105e-21),3) == 7.37

print('ALL OK')
