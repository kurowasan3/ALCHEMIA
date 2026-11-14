# 検証スクリプト ORG-DET-033  python3 -I この.py
ITEM = '/home/claude/chemx/todai/outh/items/ORG-DET-033/item.json'
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

from math import sqrt, ceil
import scipy.stats as st
def nmin(a1, a2):
    r = 2*(sqrt(a1*(1-a1)) + sqrt(a2*(1-a2)))/(a1 - a2)
    return r*r
assert ceil(nmin(0.5, 0.4)) == 392
res = {pH: nmin(0.5, a) for pH, a in ((4, 0.476), (5, 0.333), (6, 0.0833), (7, 0.0098))}
print(res)
assert sig(res[4], 2) == 6900 and ceil(res[5]) == 136 and ceil(res[6]) == 14 and ceil(res[7]) == 6
for pH, a in ((4, 0.476), (5, 0.333), (6, 0.0833), (7, 0.0098)):
    assert sig(1/(2 + 1e-5/10**-pH), 3 if pH < 7 else 2) == a
n = 392; c = 0.45*n
print('重なり', st.binom.cdf(c, n, 0.5), 1 - st.binom.cdf(c, n, 0.4))
assert st.binom.cdf(c, n, 0.5) < 0.03
sP, sQ = sqrt(60*.25), sqrt(60*.2*.8)
assert sig(2*(sP + sQ), 3) == 13.9 and 18 >= 2*(sP + sQ)

print('ALL OK')
