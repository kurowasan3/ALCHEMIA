# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem4/out/items/THE-MOL-035/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        if ' + ' not in b and not re.search(r'(^|\s)\d+\s*[A-Z(\[]', b): continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

import sympy as sp
def solve(species, sides):
    # species: list of (dict元素:数, 電荷)
    els = sorted({e for d, q in species for e in d})
    rows = [[d.get(e, 0) * s for (d, q), s in zip(species, sides)] for e in els]
    rows.append([q * s for (d, q), s in zip(species, sides)])
    ns = sp.Matrix(rows).nullspace()
    v = ns[0]; v = v / min(abs(t) for t in v if t != 0)
    l = sp.ilcm(*[sp.fraction(t)[1] for t in v]); v = v * l
    if v[0] < 0: v = -v
    return len(ns), list(v)
eq1 = [({'Zn':1},0), ({'O':1,'H':1},-1), ({'H':2,'O':1},0), ({'Zn':1,'O':4,'H':4},-2), ({'H':2},0)]
eq2 = [({'Cu':1},0), ({'H':1},1), ({'N':1,'O':3},-1), ({'Cu':1},2), ({'N':1,'O':1},0), ({'H':2,'O':1},0)]
eq3 = [({'Mn':1,'O':4},-1), ({'H':1},1), ({'Fe':1},2), ({'Mn':1},2), ({'Fe':1},3), ({'H':2,'O':1},0)]
r1 = solve(eq1, [1,1,1,-1,-1]); r2 = solve(eq2, [1,1,1,-1,-1,-1]); r3 = solve(eq3, [1,1,1,-1,-1,-1])
print(r1, r2, r3)
assert r1 == (1, [1,2,2,1,1]) and r2 == (1, [3,8,2,3,2,4]) and r3 == (1, [1,8,5,1,5,4])
# 電荷の条件を除くと (2) は一意に決まらない
els = ['Cu','H','N','O']
rows = [[d.get(e,0)*s for (d,q),s in zip(eq2,[1,1,1,-1,-1,-1])] for e in els]
assert len(sp.Matrix(rows).nullspace()) == 2
print('電荷の条件なしでは零空間 2 次元 → OK')
assert check_equation('Fe^2+ + Cl2 -> Fe^3+ + 2Cl-') is not None
assert check_equation('2Fe^2+ + Cl2 -> 2Fe^3+ + 2Cl-') is None
print('OK')
