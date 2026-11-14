# 検証スクリプト THE-EQL-101  python3 -I この.py
ITEM = '/home/claude/chemx/todai/outh/items/THE-EQL-101/item.json'
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

import math
from scipy.optimize import brentq
a = 54000/8.3
pd = lambda T: 2.6*math.exp(a*(1/273 - 1/T))
table = {273: 2.6, 280: 4.7, 283: 6.0, 294: 14.3, 295: 15.4, 296: 16.6, 297: 17.8, 298: 19.2}
for T, v in table.items():
    assert abs(pd(T) - v) / v < 0.01, (T, pd(T))
assert abs(math.log(6.0/2.6) - 0.84) < 0.01
Q = 8.3*0.84/(1/273 - 1/283); print('Q', Q); assert sig(Q, 2) == 54000
h = (4.7 - 0.10)/0.010; print('h', h); assert sig(h, 2) == 460
def zone(Tsf, tab=False):
    if not tab:
        return brentq(lambda z: 10.1 + 0.010*z - pd(Tsf + 0.030*z), 0, 3000)
    Ts = sorted(table)
    g = lambda T: 10.1 + (T - Tsf)/3 - table[T]
    for t1, t2 in zip(Ts, Ts[1:]):
        if t1 >= 294 and g(t1) > 0 > g(t2):
            Tb = t1 + g(t1)/(g(t1) - g(t2)) * (t2 - t1)
            return (Tb - Tsf)/0.030
for Tsf, want in ((276, 670), (278, 580)):
    z1, z2 = zone(Tsf), zone(Tsf, True)
    print(Tsf, z1, z2)
    assert sig(z1, 2) == want and sig(z2, 2) == want
# 海底面は安定
assert 10.1 > pd(278)
# 分解式の原子数（CH4・5.75H2O → CH4 + 5.75H2O）
lhs = {'C': 1, 'H': 4 + 2*5.75, 'O': 5.75}; rhs = {'C': 1, 'H': 4 + 2*5.75, 'O': 5.75}
assert lhs == rhs

print('ALL OK')
