# 検証スクリプト THE-SOL-044  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-SOL-044/item.json'
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

from sympy import symbols, solve, Rational as Q, nsimplify, limit, oo
b, n = symbols('b n')
s = solve([n - Q(3,10)*(90-b), n - Q(2,10)*(117-b)], [b, n]); print(s)
bv, nv = float(s[b]), float(s[n]); assert abs(bv-36.0)<1e-9 and abs(nv-16.2)<1e-9
assert sig(nv*1e-15,3) == 1.62e-14
c = nv/(150-bv); nacl = c/2; g = nacl*0.100*58.5; print('溶血', c, nacl, g)
assert sig(nacl,2) == 0.071 and sig(g,2) == 0.42
glc = 0.300*0.100*180; assert sig(glc,2) == 5.4
V = lambda cout: bv + nv/cout
assert abs(V(0.300)-90.0) < 1e-9
v1 = V(0.400); v2 = V(0.200); print(v1, v2); assert sig(v1,3) == 76.5 and sig(v2,3) == 117
c_ = symbols('c', positive=True); assert limit(bv + nv/c_, c_, 0, '+') == oo   # 尿素のみ → 溶血
for k in ('36.0','1.62\\times10^{-14}','7.1\\times10^{-2}','0.42','5.4','76.5','117'): assert k in it['answer'], k

print('全検証OK:', ITEM)
