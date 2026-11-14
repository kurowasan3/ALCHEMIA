# 検証スクリプト THE-BND-062  python3 -I この.py
ITEM = '/home/claude/chemx/todai/out/items/THE-BND-062/item.json'
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
ev = 1.602
raw = {'H': (13.60, 0.754), 'Na': (5.14, 0.548), 'F': (17.42, 3.40), 'Cl': (12.97, 3.61)}
tab = {'H': (21.8, 1.2), 'Na': (8.2, 0.9), 'F': (27.9, 5.4), 'Cl': (20.8, 5.8)}
for k, (i, e) in raw.items():
    assert abs(i*ev - tab[k][0]) < 0.1 and abs(e*ev - tab[k][1]) < 0.1, k
chi = {k: (i + e)/2 for k, (i, e) in tab.items()}
print(chi)
assert sig(chi['Cl'], 3) == 13.3
assert sorted(chi, key=chi.get, reverse=True) == ['F', 'Cl', 'H', 'Na']
assert tab['Cl'][1] > tab['F'][1] and chi['F'] > chi['Cl']
IX, EX, IZ, EZ = sp.symbols('I_X E_X I_Z E_Z')
assert sp.simplify((IX - EZ) - (IZ - EX) - ((IX + EX) - (IZ + EZ))) == 0
need = tab['Na'][0] - tab['Cl'][1]; assert sig(need, 2) == 2.4
k = 8.988e9*(1.602e-19)**2; assert sig(k, 2) == 2.3e-28
assert sig(2.3e-28/2.4e-19, 2) == 9.6e-10
d = 3.6e-30/1.27e-10/1.60e-19; print(d); assert sig(d, 2) == 0.18

print('ALL OK')
