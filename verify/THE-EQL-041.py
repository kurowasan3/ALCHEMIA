# 検証スクリプト（g04）  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-EQL-041/item.json'
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
Kp, ps = 640, 40
conv = 8.3/(8.3e3*300)        # mol/Pa
assert abs(3.0e-3/conv - 900) < 1e-6
def nosolid(H, I):
    x = symbols('x')
    sols = [s for s in solve(Eq(4*x**2, Kp*(H-x)*(I-x)), x) if 0 < s < min(H, I)]
    assert len(sols) == 1
    x = sols[0]; return H-x, I-x, 2*x
def withsolid(H):
    y = symbols('y', positive=True)
    y = solve(Eq(y**2, Kp*ps*(H - y/2)), y)[0]
    return H - y/2, ps, y
# (1)
h, i, hi = nosolid(900, 1800); print('固体なし仮定 p(I2)=', float(i)); assert i > ps
h, i, hi = withsolid(900)
chk('(1) H2', h, 100); chk('(1) I2', i, 40); chk('(1) HI', hi, 1600)
gasI = i + hi/2; chk('(1) solid mol', (1800 - gasI)*conv, 3.2e-3)
# (2)
assert gasI*conv > 1.5e-3   # 固体ありの仮定は矛盾
h, i, hi = nosolid(900, 450); assert i < ps
chk('(2) I2', i, 2.8); print('近似値', 900**2/(Kp*450))
chk('(2) approx', 900**2/(Kp*450), 2.8)
# (3)
chk('(3) min I2', gasI*conv, 2.8e-3)
# (4b) 水素を 1.0e-4 mol 加えた例
h2, _, hi2 = withsolid(900 + 1.0e-4/conv); assert hi2 > hi and h2 > 100; print('OK (4b) HI', float(hi2))
print('THE-EQL-041 すべてOK')
