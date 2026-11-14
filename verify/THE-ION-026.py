# 検証スクリプト（g04）  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-ION-026/item.json'
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
import math
Ka, Kw = 1.0e-6, 1.0e-14
def exact_h(C):
    f = lambda h: h - C*Ka/(Ka+h) - Kw/h
    lo, hi = 1e-12, 1.0
    for _ in range(300):
        mid = (lo*hi)**0.5
        if f(mid) > 0: hi = mid
        else: lo = mid
    return mid
h1 = exact_h(1.0e-2); print('(1) pH', -math.log10(h1)); assert round(-math.log10(h1), 1) == 4.0
h2 = exact_h(6.0e-6); print('(2) h', h2); assert round(-math.log10(h2), 2) == 5.70
chk('(2) ratio', math.sqrt(6.0e-12)/h2, 1.2)
h = 2.0e-7; C = (h - Kw/h)*(Ka+h)/Ka; chk('(3) C', C, 1.8e-7)
h3 = exact_h(1.8e-7); print('(3) 逆に pH', -math.log10(h3)); assert round(-math.log10(h3), 2) == 6.70
chk('(3) alpha', Ka/(Ka+h), 0.83)
Ca = h + h*h/Ka; Cb = h - Kw/h
chk('(4a)', Ca, 2.4e-7); chk('(4b)', Cb, 1.5e-7)
assert abs(Ca/C - 1) > abs(Cb/C - 1)
print('THE-ION-026 すべてOK')
