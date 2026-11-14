ITEM = '/home/claude/chemx/tokodai/out/items/THE-BND-078/item.json'
# 検証スクリプト（東工大 g04）  python3 -I この.py
import sys, json, re, math
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
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))

def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')

nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na, ('小問数', nq, na)
print('小問数OK', nq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
# ---- 個別の検証 ----

import numpy as np, itertools
from collections import Counter
def shells(basis, n=3):
    pts = []
    for i,j,k in itertools.product(range(-n,n+1), repeat=3):
        for b in basis:
            pts.append(np.array([i,j,k])+np.array(b))
    d = [round(float(np.linalg.norm(p)),6) for p in pts if np.linalg.norm(p) > 1e-9]
    c = Counter(d); ks = sorted(c)
    return [(k,c[k]) for k in ks[:4]]
fcc = shells([(0,0,0),(.5,.5,0),(.5,0,.5),(0,.5,.5)])
bcc = shells([(0,0,0),(.5,.5,.5)])
print(fcc, bcc)
assert fcc[0][1] == 12 and abs(fcc[0][0]-2**.5/2) < 1e-6 and fcc[1] == (1.0, 6)
assert bcc[0][1] == 8 and abs(bcc[0][0]-3**.5/2) < 1e-6 and bcc[1] == (1.0, 6)
chk('fcc比', fcc[1][0]/fcc[0][0], 1.41, 3)
chk('bcc比', bcc[1][0]/bcc[0][0], 1.15, 3)

print('ALL OK')
