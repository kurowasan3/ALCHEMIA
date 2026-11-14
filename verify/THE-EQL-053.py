# 検証スクリプト THE-EQL-053  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-EQL-053/item.json'
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

from sympy import Rational as Q, symbols, Poly, diff, log, simplify, solve
import random
K = Q(1,10)**2/(Q(6,10)*Q(3,10)**3); print('K', K, float(K)); assert sig(float(K),3) == 0.617
q = lambda nN2, nH2, nNH3: Q(nNH3)**2*(nN2+nH2+nNH3)**2/(Q(nN2)*Q(nH2)**3)
assert q(6,3,1) == K
q2 = q(7,3,1); print('直後', q2); assert q2 == Q(121,189) and q2 > K and sig(float(q2),3) == 0.640
x = symbols('x')
f = (1-2*x)**2*(11+2*x)**2 - K*(7+x)*(3+3*x)**3
roots = [r for r in Poly(f, x).nroots() if r.is_real and 0 < r < 0.5]; print(roots)
assert len(roots) == 1; nh3 = 1-2*float(roots[0]); assert sig(nh3,3) == 0.989
# 体積一定: Kc の反応商は n_N2 増で減少
V = symbols('V', positive=True); a,b,c = symbols('a b c', positive=True)
qc = c**2*V**2/(a*b**3); assert simplify(diff(log(qc), a)) == -1/a
# 全圧一定の微分
lq = 2*log(c) + 2*log(a+b+c) - log(a) - 3*log(b)
dN2 = simplify(diff(lq, a)); dH2 = simplify(diff(lq, b)); print(dN2, dH2)
assert simplify(dN2 - (2/(a+b+c) - 1/a)) == 0 and simplify(dH2 - (2/(a+b+c) - 3/b)) == 0
for _ in range(1000):
    nN, nH, nA = [random.uniform(0.01, 10) for _ in range(3)]; N = nN+nH+nA
    assert 2/N - 3/nH < 0
    D = random.uniform(0.001, 20)
    ratio = ((N+D)/N)**2 * nN/(nN+D)
    assert (ratio > 1) == ((nN/N)*(2 + D/N) > 1)
t1 = 1/0.45 - 2; t2 = 1/0.25 - 2; print(t1, t2); assert sig(t1,3) == 0.222 and sig(t2,3) == 2.00
assert check_equation('N2 + 3H2 -> 2NH3') is None

print('全検証OK:', ITEM)
