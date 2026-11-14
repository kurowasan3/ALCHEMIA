# 検証スクリプト THE-GAS-071   python3 -I この.py
ITEM = '/home/claude/chemk/out/items/THE-GAS-071/item.json'
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
NA, R = 6.02e23, 8.31
d = lambda P, T: P/(R*T)*NA/1e6
d0 = d(1.013e5, 273); print('(1)', d0); assert sig(d0,3) == 2.69e19
c0 = d0*1.90e-6; print('(2)', c0); assert sig(c0,3) == 5.11e13
d1 = d(5.00e4, 253); c1 = d1*1.90e-6; print('(3)', c1, c1/c0); assert sig(c1,3) == 2.72e13 and sig(c1/c0,2) == 0.53
m = 1.013e5/(R*298)*40.0e-9*48.0*1e6; print('(4)', m); assert sig(m,3) == 78.5

print('検証OK')
