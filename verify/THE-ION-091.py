# 検証スクリプト THE-ION-091  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/out/items/THE-ION-091/item.json'
import sys, json, re
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    if n_items: assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()


import sympy as sp
KP, KQ = 1.0e-4, 1.0e-6
a = lambda K, pH: K/(K+10**-pH)
assert abs(a(KP,4.0)-0.5) < 1e-12
assert sig(a(KP,5.0),2) == 0.91 and sig(a(KQ,5.0),2) == 0.091
d = {pH: a(KP,pH)-a(KQ,pH) for pH in (3.0,5.0,7.0,9.0)}
print(d)
assert max(d, key=d.get) == 5.0
assert [round(d[p],3) for p in (3.0,5.0,7.0,9.0)] == [0.090, 0.818, 0.090, 0.001]
h = sp.symbols('h', positive=True)
f = KP/(KP+h) - KQ/(KQ+h)
sol = [s for s in sp.solve(sp.diff(f,h), h) if s.is_positive]
print('最大の [H+]', sol); assert abs(-sp.log(sol[0],10) - 5.0) < 1e-9
print('OK')
