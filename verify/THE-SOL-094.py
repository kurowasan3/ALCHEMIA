# 検証スクリプト（東工大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/outh/items/THE-SOL-094/item.json'
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

from sympy import symbols, solve, nsolve
V, p, q = symbols('V p q', positive=True)
o = 0.049*2.0*10; n = 0.024*8.0*10
assert sig(o+n,3) == 2.90 and sig(o/(o+n)*100,3) == 33.8
r = 0.049/0.024; assert sig(r,3) == 2.04
assert 0.25*r**3 < 4 <= 0.25*r**4
s = solve([p*(V+0.49)-1, q*(V+0.24)-4, p+q-10], [V,p,q], dict=True); assert len(s) == 1; s = s[0]
assert sig(float(s[V]),3) == 0.225 and sig(float(s[p]),3) == 1.40
dO, dN = 0.49*float(s[p]), 0.24*float(s[q]); assert sig(dO/(dO+dN)*100,3) == 24.9
assert sig(float(s[p])/10*100,3) == 14.0

print('ALL OK')
