# 検証スクリプト（g04）  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-EQL-040/item.json'
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

# 全圧と A 分圧の関係
pA, pA0, pB0 = symbols('pA pA0 pB0', positive=True)
P = pA + pB0 + (pA0 - pA)/2
pA_of = lambda Pt, a0, b0: solve(Eq(P.subs({pA0: a0, pB0: b0}), Pt), pA)[0]
# (1)
a, b = symbols('a b', positive=True)
s = solve([Eq(a + b, 2.0e4), Eq(b + a/2, 1.5e4)], [a, b])
chk('(1) pA0', s[a], 1.0e4); chk('(1) pB0', s[b], 1.0e4)
# (2)
p1 = pA_of(3.0e4, 4.0e4, 0); assert abs(p1 - 2.0e4) < 1e-6
p2 = pA_of(1.75e4, s[a], s[b]); assert abs(p2 - 0.5e4) < 1e-6
t1, t2 = 50, 200
# n=1 なら半減期は等しいはず → 矛盾
assert t1 != t2
# n=2: t1/2 = 1/(k p0)
k1 = 1/(t1*4.0e4); k2 = 1/(t2*float(s[a])); assert abs(k1-k2) < 1e-15
chk('(3) k', k1, 5.0e-7)
# (4)
pa = 1/(1/4.0e4 + k1*150); chk('(4) P', (4.0e4 + pa)/2, 2.5e4)
# 数値積分でも確認
import math
def integrate(p0, k, T, n=2, dt=0.01):
    p = p0; t = 0
    while t < T - 1e-9:
        f = lambda q: -k*q**n
        k1_ = f(p); k2_ = f(p+dt/2*k1_); k3_ = f(p+dt/2*k2_); k4_ = f(p+dt*k3_)
        p += dt/6*(k1_+2*k2_+2*k3_+k4_); t += dt
    return p
assert abs(integrate(4.0e4, k1, 150) - 1.0e4) < 1
# (5)
pa5 = pA_of(1.60e4, s[a], s[b]); chk('(5) pA', pa5, 2.0e3)
t5 = (1/pa5 - 1/s[a])/k1; chk('(5) t', t5, 8.0e2)
print('THE-EQL-040 すべてOK')
