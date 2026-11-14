ITEM = '/home/claude/chemx/tokodai/out/items/THE-THC-056/item.json'
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

import itertools
q = {'H2':286,'CH4':891,'C3H8':2219,'C4H10':2878}
c = {'H2':0,'CH4':1,'C3H8':3,'C4H10':4}
chk('(1)CH4', 1/891*1000, 1.12, 3)
chk('(1)C4H10', 4/2878*1000, 1.39, 3)
res = {}
for a,b in itertools.combinations(q,2):
    lo,hi = sorted([q[a],q[b]])
    if not (lo <= 1200 <= hi): res[(a,b)] = None; continue
    if q[a] > q[b]: a,b = b,a
    y = (1200-q[a])/(q[b]-q[a])
    res[(a,b)] = (1-y)*c[a]+y*c[b]
print(res)
assert res[('H2','CH4')] is None and res[('C3H8','C4H10')] is None
best = min((v,k) for k,v in res.items() if v is not None)
assert best[1] == ('H2','C4H10')
chk('(3)', best[0], 1.41, 3)
# 3成分以上も含めた線形計画
from scipy.optimize import linprog
r = linprog([0,1,3,4], A_eq=[[1,1,1,1],[286,891,2219,2878]], b_eq=[1,1200], bounds=[(0,None)]*4)
print('LP', r.fun, r.x)
chk('LP', r.fun, 1.41, 3)
assert abs(r.fun-best[0]) < 1e-9
for e in ['2H2 + O2 -> 2H2O','CH4 + 2O2 -> CO2 + 2H2O','C3H8 + 5O2 -> 3CO2 + 4H2O','2C4H10 + 13O2 -> 8CO2 + 10H2O']:
    assert check_equation(e) is None

print('ALL OK')
