# 検証スクリプト THE-ATM-032  python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-ATM-032/item.json'
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
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|\$', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n += 1
    print('反応式の検査', n, '本 OK')
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

from fractions import Fraction as Fr
from sympy import symbols, simplify, pi
n, h, k, m, e, Z = symbols('n h k m e Z', positive=True)
rH = n**2*h**2/(4*pi**2*k*m*e**2); EH = -2*pi**2*k**2*m*e**4/(n**2*h**2)
rZ = rH.subs(e**2, Z*e**2); EZ = -2*pi**2*k**2*m*(Z*e**2)**2/(n**2*h**2)
assert simplify(rZ - rH/Z) == 0 and simplify(EZ - EH*Z**2) == 0
assert sig(0.0529/2, 3) == 0.0265 and abs(13.6*4 - 54.4) < 1e-9
# (2) 全列挙
def lev(Zv, a, b): return Zv*Zv*(Fr(1, a*a) - Fr(1, b*b))
H = {}
for a in range(1, 41):
    for b in range(a+1, 41): H.setdefault(lev(1, a, b), []).append((b, a))
coin = []
for a in range(1, 21):
    for b in range(a+1, 21):
        E = lev(2, a, b)
        if E in H: coin.append(((b, a), H[E]))
for a in range(1, 10):
    for b in range(a+1, 11):
        assert lev(2, 2*a, 2*b) in H
odd = [c for c in coin if not (c[0][0] % 2 == 0 and c[0][1] % 2 == 0)]
print('偶然の一致', odd); assert ((9, 6), [(36, 4)]) in odd
assert not any(Fr(1, 4) - Fr(7, 36) == Fr(1, b*b) for b in range(1, 1000))
assert ((6, 4), [(3, 2)]) in coin
E43 = lev(2, 3, 4); assert E43 not in H
lam = 1240/(13.6*float(E43)); print('He+ 4->3', 13.6*float(E43), lam)
assert sig(lam, 3) == 469
# (3)
hits_lim = []; hits = []
for Zv in range(1, 11):
    for nL in range(1, 11):
        lim = 1240/(13.6*Zv*Zv/nL**2)
        first = 1240/(13.6*float(lev(Zv, nL, nL+1)))
        if abs(lim/91.2 - 1) < 0.005:
            hits_lim.append((Zv, nL))
            if abs(first/208 - 1) < 0.005: hits.append((Zv, nL))
print('limit only', hits_lim, 'both', hits)
assert all(a == b for a, b in hits_lim) and len(hits_lim) > 1 and hits == [(3, 3)]
for nL, v in [(1, 0.750), (2, 0.556), (3, 0.4375), (4, 0.360)]:
    assert abs(float(1 - Fr(nL*nL, (nL+1)**2)) - v) < 6e-4
assert sig(91.2/208, 3) == 0.438
# (4)
assert abs(-13.6*9/9 + 13.6) < 1e-12 and sig(0.0529*9/3, 3) == 0.159

print('ALL OK')
