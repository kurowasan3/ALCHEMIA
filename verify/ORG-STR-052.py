# 検証スクリプト ORG-STR-052  python3 -I この.py
ITEM = '/home/claude/chemx/todai/out/items/ORG-STR-052/item.json'
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

import numpy as np
b, r, th = 1.54, 1.43, np.radians(109.5)
table = {0: 2.5, 60: 2.8, 120: 3.4, 180: 3.7}
for phi, v in table.items():
    p = np.radians(phi)
    O1 = np.array([-r*np.cos(np.pi-th), r*np.sin(np.pi-th), 0])
    O2 = np.array([b + r*np.cos(np.pi-th), r*np.sin(np.pi-th)*np.cos(p), r*np.sin(np.pi-th)*np.sin(p)])
    d = np.linalg.norm(O2-O1)
    assert round(d, 1) == v, (phi, d)
    print(phi, round(d, 2))
hb = [k for k, v in table.items() if v <= 3.0]
assert hb == [0, 60]

print('ALL OK')
