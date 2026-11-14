# 検証スクリプト（g04）  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-ION-025/item.json'
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
Kw = 1e-14
def exact_h(Ka, C0, V0, cb, v):
    V = V0 + v; Na = cb*v/V; C = C0*V0/V
    f = lambda h: h + Na - C*Ka/(Ka + h) - Kw/h
    lo, hi = 1e-13, 1.0
    for _ in range(200):
        mid = (lo*hi)**0.5
        if f(mid) > 0: hi = mid
        else: lo = mid
    return mid
# 与えたデータが真のモデルと整合
chk('data v=4', exact_h(2.0e-5, 0.080, 20.0, 0.100, 4.00), 6.0e-5)
chk('data v=12', exact_h(2.0e-5, 0.080, 20.0, 0.100, 12.00), 6.7e-6)
Ve = symbols('Ve', positive=True)
Ves = solve(Eq(6.0e-5/6.7e-6, (Ve-4)/4/((Ve-12)/12)), Ve)[0]
chk('(2) Ve', Ves, 16.0, 3)
Ka = 6.0e-5*4/(Ves-4); chk('(2) Ka', Ka, 2.0e-5)
chk('(2) C', 0.100*Ves/20.0, 8.0e-2)
# (3)
Ka = 2.0e-5; Vev = 16.0
chk('(3) approx', Ka*(Vev-0.40)/0.40, 7.8e-4)
h = symbols('h', positive=True)
Na = 0.100*0.40/20.4; Ct = 0.100*Vev/20.4
hs = [s for s in solve(Eq(h*(Na+h), Ka*(Ct-Na-h)), h)][0]
chk('(3) exact', hs, 5.9e-4)
chk('(3) exact full', exact_h(2.0e-5, 0.080, 20.0, 0.100, 0.40), 5.9e-4)
print('THE-ION-025 すべてOK')
