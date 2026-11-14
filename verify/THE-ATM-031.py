ITEM = '/home/claude/chemk/outh/items/THE-ATM-031/item.json'
# 検証スクリプト（グループ g04）  python3 -I この.py
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

import sympy as sp
k, b, m, h, c = sp.symbols('k b m h c', positive=True)
L = (2*k+4)*b
E = lambda n: n**2*h**2/(8*m*L**2)
dE = sp.simplify(E(k+3) - E(k+2))
assert sp.simplify(dE - (2*k+5)*h**2/(8*m*L**2)) == 0
lam = sp.simplify(h*c/dE)
assert sp.simplify(lam - 8*m*c*(2*k+4)**2*b**2/((2*k+5)*h)) == 0
vals = {m: 9.11e-31, h: 6.63e-34, c: 3.00e8}
bb = sp.sqrt(9*6.63e-34*416e-9/(8*9.11e-31*3.00e8*64))
print('b =', float(bb)); assert sig(float(bb)*1e9, 3) == 0.133
lamf = sp.lambdify(k, lam.subs(vals).subs(b, 0.133e-9))
ratio = lambda kk: (2*kk+4)**2/(2*kk+5) / (64/9) * 416
for kk, exp in [(1, 3.0e2), (3, 5.3e2), (4, 6.5e2)]:
    a1, a2 = lamf(kk)*1e9, ratio(kk)
    print(kk, a1, a2); assert sig(a1, 2) == exp and sig(a2, 2) == exp
def color(l):
    if l < 400: return '無色'
    for lo, hi, col in [(400,435,'黄緑'),(435,480,'黄'),(480,500,'橙〜赤'),(500,560,'赤紫'),(560,610,'青'),(610,750,'青緑')]:
        if lo <= l < hi: return col
assert [color(ratio(kk)) for kk in (1,3,4)] == ['無色','赤紫','青緑']
assert [color(lamf(kk)*1e9) for kk in (1,3,4)] == ['無色','赤紫','青緑']
inc = 16*9.11e-31*3.00e8*(0.133e-9)**2/6.63e-34
print('増分', inc); assert sig(inc*1e9, 2) == 1.2e2
first = [kk for kk in range(1, 9) if ratio(kk) > 700][0]; first2 = [kk for kk in range(1, 9) if lamf(kk)*1e9 > 700][0]
print('700nm超', first, first2, ratio(5), lamf(5)*1e9); assert first == first2 == 5
print('差', [round(ratio(kk+1)-ratio(kk),1) for kk in range(1,8)])

print('全検証OK:', ITEM)
