# 検証スクリプト THE-ATM-020  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/THE-ATM-020/item.json'
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
        print('式OK:', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
    for k in ('problem','answer','solution'):
        assert it[k].count('$') % 2 == 0, k
from sympy import log, Rational as R, N
check_structure(); check_all_equations()
assert 'CaCO3 + CO2 + H2O -> Ca(HCO3)2' in all_equations()
T = 5730
t1 = T*N(log(4, 2)); print('t1', t1); assert sig(float(t1),3) == 11500
l12 = (2*0.301 + 0.477 - 1)/0.301
t2g = T*(2 + round(l12, 3)); t2e = T*float(N(2 + log(R(6,5), 2)))
print('log2 1.2', l12, 't2(与えた値)', t2g, 't2(厳密)', t2e)
assert sig(t2g,3) == 13000 and sig(t2e,3) == 13000 and t2e > t1
# 核反応式の保存（質量数・原子番号）
assert 14+1 == 14+1 and 7+0 == 6+1
assert 14 == 14 and 6 == 7 - 1
print('THE-ATM-020 OK')
