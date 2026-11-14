# 検証スクリプト（g04）  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-ACB-018/item.json'
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
def comp(V1, V2):
    a, b, c = symbols('a b c', nonnegative=True)
    res = []
    for zero in (a, c):
        s = solve([Eq(0.100*V1, a + b), Eq(0.100*V2, b + c), Eq(zero, 0)], [a, b, c], dict=True)
        for d in s:
            if all(v >= -1e-12 for v in d.values()): res.append(tuple(float(d[k]) for k in (a, b, c)))
    res = sorted(set(res)); assert len(res) == 1, res   # 一意
    return res[0]
X = comp(6.0, 14.0); Y = comp(15.0, 5.0); Z = comp(0, 8.0)
print('X', X, 'Y', Y, 'Z', Z)
chk('X Na2CO3', X[1]*1e-3, 6.0e-4); chk('X NaHCO3', X[2]*1e-3, 8.0e-4); assert X[0] == 0
chk('Y NaOH', Y[0]*1e-3, 1.0e-3); chk('Y Na2CO3', Y[1]*1e-3, 5.0e-4); assert Y[2] == 0
chk('Z NaHCO3', Z[2]*1e-3, 8.0e-4); assert Z[0] == 0 and Z[1] == 0
Na = Y[0] + 2*Y[1]; chk('(3) c', Na*1e-3/0.0100, 0.20); chk('(3) CO2', Y[1]*1e-3, 5.0e-4)
b, c = symbols('b c', positive=True)
s = solve([Eq(c, 2*b), Eq(2*b + c, Na)], [b, c])
chk('(4) V1', s[b]/0.100, 5.0); chk('(4) V2', (s[b]+s[c])/0.100, 15.0, 3); chk('(4) CO2', (s[b]+s[c])*1e-3, 1.5e-3)
assert abs((s[b]+s[c]+s[b])/0.100 - 20.0) < 1e-9
print('THE-ACB-018 すべてOK')
