# 検証スクリプト THE-EQL-033  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-EQL-033/item.json'
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
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b:
                eqs.add(b)
    return sorted(eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
from sympy import Rational as R, sqrt, nsimplify
check_structure()
eqs = all_equations(); print(eqs)
for e in eqs: assert check_equation(e) is None, e
assert '(1)' not in re.sub(r'^\(1\)', '', it['answer'])
tab = [(R(4), R(60,10000)), (R(16), R(12,1000)), (R(36), R(18,1000)), (R(64), R(24,1000))]
ka = {x/p for p, x in tab}; kb = {x**2/p for p, x in tab}
print('x/p', [float(v) for v in ka], 'x^2/p', kb)
assert len(ka) > 1 and kb == {R(9, 10**6)}
Kb = R(9, 10**6)
x = lambda p: sqrt(Kb*p)
x2 = x(25); assert x2 == R(15,1000)
ratio = (x(1600) - x2)/(x(400) - x2); print('ratio', float(ratio)); assert sig(float(ratio), 3) == 2.33
lin = R(1600-25, 400-25); assert sig(float(lin), 3) == 4.20
pH2 = 400*R(1,4); assert pH2 == 100
r3 = (x(pH2) - x2)/(x(400) - x2); assert sig(float(r3), 2) == 0.33
for p, simple, corr in ((100, 0.0300, 0.0291), (10**5, 0.949, 0.487)):
    xs = float(x(p)); y = xs; xc = y/(1+y)
    # 補正式を直接解いて確認
    from sympy import symbols, solve
    X = symbols('X', positive=True)
    sol = [s for s in solve((X/(1-X))**2 - Kb*p, X) if 0 < s < 1]
    print(p, xs, xc, [float(s) for s in sol])
    assert sig(xs, 3) == simple and sig(xc, 3) == corr and abs(float(sol[0]) - xc) < 1e-12
print('x=1 となる圧力', float(1/Kb))
print('全検証OK:', ITEM)
