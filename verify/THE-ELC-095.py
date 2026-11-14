# 検証スクリプト THE-ELC-095   python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-ELC-095/item.json'
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
        print('式', e, '→', 'OK' if r is None else r)
        assert r is None, e
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()
# ---- 個別の検証 ----
from math import log10
s = 0.0592/2
ECu = lambda c: 0.337 + s*log10(c)
ECo = lambda c: -0.277 + s*log10(c)
ENi = lambda c: -0.257 + s*log10(c)
P0, Q0 = ECu(1)-ECo(0.1), ECu(1)-ENi(1e-3)
print('(1)', P0, Q0); assert (sig(P0,3), sig(Q0,3)) == (0.644, 0.683)
import sympy as sp
x = sp.symbols('x', positive=True)
r = 10**(-0.020/s); xs = sp.solve(sp.Eq((1e-3+x), r*(0.1+x)), x)[0]
t = float(2*xs*9.65e4/0.5/3600); print('(2) x', xs, 't/h', t); assert sig(t,3) == 2.73
P1, Q1 = ECu(0.95)-ECo(0.15), ECu(0.95)-ENi(0.051)
print('(3)', P1, Q1); assert (sig(P1,3), sig(Q1,3)) == (0.638, 0.632)
cmax = 10**((ECu(0.95)-0.650+0.257)/s)
V = 0.05/(cmax-1e-3); print('(4)', cmax, V); assert sig(V,3) == 4.47
# V でのチェック：起電力がちょうど 0.650
assert abs(ECu(0.95)-ENi(1e-3+0.05/V) - 0.650) < 1e-9

print('検証OK')
