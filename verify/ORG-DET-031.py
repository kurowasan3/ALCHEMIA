# 検証スクリプト ORG-DET-031  python3 -I この.py
ITEM = '/home/claude/chemx/todai/out/items/ORG-DET-031/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0:
            return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

def check_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if should_check(b):
                r = check_equation(clean_ce(b))
                assert r is None, (b, r)
                n += 1
    print('反応式', n, '本 OK')

def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    if n_items:
        assert n_items == n_ans, (n_items, n_ans)
    else:
        assert n_ans == 0
    print('小問数', n_items, '答え数', n_ans)

def check_svgs():
    import os
    d = os.path.dirname(ITEM)
    for k in ('problem_fig', 'solution_fig'):
        f = it[k]
        if f:
            s = open(os.path.join(d, f), encoding='utf-8').read()
            assert 'viewBox' in s and '<script' not in s and 'foreignObject' not in s
            assert not re.search(r'\son[a-z]+=', s) and 'href="http' not in s and 'base64' not in s
    print('SVG OK')

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

check_structure(); check_equations(); check_svgs()

# ---- 個別の検証 ----

import sympy as sp
from math import comb
def simulate(a, n=4):
    water = [0.0]*(n+1); org = [0.0]*(n+1)
    tot = [0.0]*(n+1); tot[0] = 1.0
    org[0] = a*tot[0]; water[0] = (1-a)*tot[0]
    for step in range(n):
        org = [0.0] + org[:-1]          # 有機層を隣へ
        tot = [w + o for w, o in zip(water, org)]
        org = [a*t for t in tot]; water = [(1-a)*t for t in tot]
    return tot
for a in (0.70, 0.25):
    s = simulate(a); b = [comb(4, k)*a**k*(1-a)**(4-k) for k in range(5)]
    assert all(abs(x - y) < 1e-12 for x, y in zip(s, b)), (s, b)
    print(a, [round(v, 4) for v in s])
X = simulate(0.70); Y = simulate(0.25)
assert [sig(v, 2) for v in X] == [0.0081, 0.076, 0.26, 0.41, 0.24]
assert Y.index(max(Y)) == 1 and sig(max(Y), 2) == 0.42
a = sp.symbols('a', positive=True)
sol = sp.solve(sp.Eq(6*a**2*(1-a)**2, 2*4*a*(1-a)**3), a)
sol = [s for s in sol if 0 < s < 1]; print(sol); assert sig(float(sol[0]), 2) == 0.57
rx = X[3] + X[4]; ry = Y[3] + Y[4]
assert sig(100*rx, 2) == 65 and sig(100*rx/(rx+ry), 2) == 93
assert sig(100*0.70/(0.95), 2) == 74

print('ALL OK')
