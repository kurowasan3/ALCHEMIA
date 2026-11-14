# 検証スクリプト THE-SOL-069  python3 -I この.py
ITEM = '/home/claude/chemx/todai/out/items/THE-SOL-069/item.json'
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
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm', b): continue
            if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', b): continue  # 一般式（validate と同じ除外）
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n += 1
    print('反応式', n, '本 OK')
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure(); check_all_equations()
import sympy as sp
R, T, C = 8.3e3, 300, 1.00e-2
P0 = C*R*T; assert sig(P0, 3) == 2.49e4
i = 2.74e4/P0; a = i - 1; assert sig(a, 2) == 0.10
Ka = C*0.10**2/0.90; assert sig(Ka, 2) == 1.1e-4
i3 = 6.97e4/P0; k = 1 + (i3-1)/0.90; assert round(k) == 3 and abs(k-3) < 0.01
assert abs(0.953/1.00e-2 - (24.3+2*35.5)) < 0.05
assert abs(40.1+71.0 - 95.3) > 10
x = sp.symbols('x', positive=True)
ad = sp.nsolve(1.00e-3*x**2/(1-x) - Ka, x, 0.3); Pd = 1.00e-3*R*T*(1+ad)
print('希釈後 α', ad, 'Π', Pd); assert Pd > 2.74e4/10 and sig(float(Pd), 2) == 3.2e3
print('ALL OK')
