ITEM = '/home/claude/chemx/tokodai/outh/items/THE-BND-084/item.json'
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
def pts_hcp(n=4):
    a1=np.array([1,0,0]); a2=np.array([0.5,3**.5/2,0]); c=np.array([0,0,(8/3)**.5])
    return np.array([(i+b[0])*a1+(j+b[1])*a2+(k+b[2])*c for i,j,k in itertools.product(range(-n,n+1),repeat=3) for b in [(0,0,0),(1/3,1/3,0.5)]])
def pts_fcc(n=4):
    return np.array([(np.array([i,j,k])+np.array(b))*2**.5 for i,j,k in itertools.product(range(-n,n+1),repeat=3) for b in [(0,0,0),(.5,.5,0),(.5,0,.5),(0,.5,.5)]])
res = {}
for nm,P in [('fcc',pts_fcc()),('hcp',pts_hcp())]:
    d2 = np.sum(P**2,axis=1); d2 = d2[d2>1e-9]
    c = Counter(np.round(d2,5)); res[nm] = c
    print(nm, sorted(c.items())[:7])
f, h = res['fcc'], res['hcp']
assert f[1.0]==12 and h[1.0]==12
assert f[2.0]==6 and h[2.0]==6 and f[3.0]==24 and h[3.0]==18
assert h[round(8/3,5)]==2 and round(8/3,5) not in f
assert h[round(11/3,5)]==12 and f[4.0]==12 and h[4.0]==6
assert f[6.0] == 8   # √6 d（体対角線方向を含む殻）
# 面心で真上方向（[111]）の最初の粒子
v = np.array([1,1,1])/3**.5
P = pts_fcc(); al = [np.linalg.norm(p) for p in P if np.linalg.norm(p)>1e-9 and np.linalg.norm(np.cross(p,v))<1e-6]
chk('(1)fcc', min(al), 2.45, 3)
chk('(1)hcp', (8/3)**.5, 1.63, 3)
within = lambda c,R2: sum(v for k,v in c.items() if k <= R2+1e-6)
assert within(f,3)==42 and within(h,3)==38
assert within(f,4)==54 and within(h,4)==56
print('連続近似(X含む):', 4/3*np.pi*8*2**.5)

print('ALL OK')
