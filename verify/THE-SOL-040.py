# 検証スクリプト THE-SOL-040  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-SOL-040/item.json'
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
def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    return re.sub(r'\{(\d*[+-])\}', r'\1', eq)
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'<=>|->|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                e = check_equation(clean_ce(b)); n += 1
                assert e is None, (b, e)
    print('反応式', n, '本 つり合いOK')
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

from sympy import Rational as Q, nsimplify, symbols, solve
RT = 8.31e3*300
pi1 = 1.27*98.0; MP = 3.00/(pi1/RT); print('M_P', MP); assert sig(MP,2) == 6.0e4
pi2 = 3.82*98.0; n = pi2/RT
MQ = 1.50/(n - 1.50/6.0e4); print('M_Q', MQ); assert sig(MQ,2) == 1.2e4
Mapp = 3.00/n; print('見かけ', Mapp); assert sig(Mapp,2) == 2.0e4 and Mapp < (6.0e4+1.2e4)/2
# 粒子数比 Q/P = 5
assert abs((1.50/1.2e4)/(1.50/6.0e4) - 5) < 1e-9
h = (pi2-196)/98.0; print('h', h); assert sig(h,3) == 1.82
assert sig(pi2,3) == 374
for k in ('6.0\\times10^{4}','1.2\\times10^{4}','2.0\\times10^{4}','1.82','374'): assert k in it['answer'], k

print('全検証OK:', ITEM)
