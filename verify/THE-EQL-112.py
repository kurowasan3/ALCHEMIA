# 検証スクリプト THE-EQL-112  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/out/items/THE-EQL-112/item.json'
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
a = sp.symbols('a', positive=True)
a1 = sp.solve(sp.Eq(3.0 - a/2, 2.4), a)[0]; print(a1); assert abs(a1-1.2) < 1e-12
n = {'CO': 2.0-a1, 'O2': 1.0-a1/2, 'CO2': a1}; P = 2.4e5; N = sum(n.values())
p = {k: v/N*P for k, v in n.items()}
Kp = p['CO2']**2/(p['CO']**2*p['O2']); print(Kp); assert sig(float(Kp),2) == 5.6e-5
dP, V, R, T = sp.symbols('dP V R T', positive=True)
assert sp.simplify(sp.solve(sp.Eq(dP, a/2*R*T/V), a)[0] - 2*dP*V/(R*T)) == 0
print('OK')
