# 検証スクリプト（g04）  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-EQL-042/item.json'
from sympy import *
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
print('反応式', len(eqs), '本を確認')

def sig(x, n=2):
    x = float(x)
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def chk(name, got, want, n=2):
    assert sig(got, n) == sig(want, n), (name, float(got), want)
    print('OK', name, sig(got, n))
x, y = Rational(75, 100), Rational(1) - Rational(75, 100)
G, M, D, A, W = 1 - x, x - y, y, Rational(150, 100) - x - y, x + y
assert A == Rational(1, 2) and G == Rational(1, 4)
K1 = M*W/(G*A); K2 = D*W/(M*A)
chk('(1) M', M, 0.50); chk('(1) D', D, 0.25); chk('(1) W', W, 1.00); chk('(1) K1', K1, 4.0); chk('(1) K2', K2, 1.0)
r = A/W; assert (G, M, D) == (G, G*K1*r, G*K1*K2*r**2); print('OK (2) 分布式')
r = 9; tot = 1 + K1*r + K1*K2*r**2
Mf, Df = K1*r/tot, K1*K2*r**2/tot; W = Mf + 2*Df; a = r*W + W
chk('(3) a', a, 19); print('a =', float(a))
# 直接連立で確認
xs, ys = symbols('xs ys', positive=True)
sol = nsolve([ (xs-ys)*(xs+ys) - K1*(1-xs)*(a-xs-ys), ys*(xs+ys) - K2*(xs-ys)*(a-xs-ys)], [xs, ys], [0.99, 0.9])
print('D/M =', float(sol[1]/(sol[0]-sol[1])))
assert abs(sol[1]/(sol[0]-sol[1]) - 9) < 1e-6
print('THE-EQL-042 すべてOK')
