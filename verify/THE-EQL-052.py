# 検証スクリプト THE-EQL-052  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-EQL-052/item.json'
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

from math import log, exp
V = lambda t: 50.0*(1 - 0.8**(t/2))
table = [0, 8.00, 14.40, 19.52, 23.62, 26.89]
for t, v in zip((2,4,6,8,10,12), table): assert abs((V(t)-V(2)) - v) < 0.006, (t, V(t)-V(2), v)
inc = [table[i+1]-table[i] for i in range(5)]; rat = [inc[i+1]/inc[i] for i in range(4)]; print(inc, rat)
assert all(abs(r-0.8) < 0.003 for r in rat)
# 二次反応の比較: 同じ最初の区間の減少率をもつ二次反応では比が増えていく
C0 = 1.0; kk = (1/0.8 - 1)/(C0*2)   # C(2)=0.8 C0 となる二次反応
C = lambda t: C0/(1 + kk*C0*t)
d2 = [C(t)-C(t+2) for t in (2,4,6,8,10)]; r2 = [d2[i+1]/d2[i] for i in range(4)]; print('二次', r2)
assert all(r2[i+1] > r2[i] for i in range(3))
k = (log(5) - 2*log(2))/2; print('k', k); assert sig(k,3) == 0.112
th = log(2)/k; assert sig(th,3) == 6.21
k_given = (1.6094 - 2*0.6931)/2; print(k_given, 0.6931/k_given); assert sig(0.6931/0.1116,3) == 6.21 and sig(k_given,3) == 0.112
vinf = 8.00/(1-0.800); assert abs(vinf-40.0) < 1e-9
n = 1.013e5*0.0400/(8.31e3*298); c0 = 2*n/0.800/0.0100; print('c0', c0); assert sig(c0,3) == 0.409
kw = log(26.89/(26.89-8.00))/2.0; print('誤った k', kw, kw/k); assert sig(kw,3) == 0.177 and kw > k
Ea = 8.31*log(2)/(1/298-1/308); print('Ea', Ea); assert sig(Ea/1000,3) == 52.9
assert sig(8.31*0.6931/(1/298-1/308)/1000,3) == 52.9
assert abs(log(1/0.64)/log(1/0.8) - 2) < 1e-12

print('全検証OK:', ITEM)
