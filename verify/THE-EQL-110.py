# 検証スクリプト（東大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/todai/outh/items/THE-EQL-110/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

# 小問の数と答えの数の一致
nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans == 0), (nprob, nans)
print('小問数', nprob, '答え', nans)

# ---- 個別の検証 ----
from sympy import symbols, nsolve, sqrt, Rational
x, s = symbols('x s', positive=True)
Kp = 3.0e10
def eq(P, s_, extra=0):
    return 27*x**4*P**2 - Kp*(1-x)*(s_-x)*(1+s_+2*x+extra)**2
x1 = nsolve(eq(1e5, 1), x, 0.6); x2 = nsolve(eq(1e6, 1), x, 0.3)
xN = nsolve(eq(1e6, 1, 2), x, 0.3); x3 = nsolve(eq(1e6, 3), x, 0.4)
print(x1, x2, xN, x3)
assert sig(float(x1), 2) == 0.63 and sig(float(x2), 2) == 0.25 and sig(float(xN), 2) == 0.32
assert abs(float(xN) - (-1 + 129**0.5)/32) < 1e-9
assert x3 > xN
xx = 0.32
Q = 27*xx**4*1e12/((1-xx)*(3-xx)*(4+2*xx)**2)
print('Q', Q); assert Q < Kp and sig(Q, 1) == 7e9
sreq = nsolve(27*0.5**4*1e12 - Kp*0.5*(s-0.5)*(s+2)**2, s, 4)
print('s', sreq); assert sig(float(sreq), 2) == 3.8
assert abs((3.8-0.5)*(3.8+2)**2 - 111.0) < 0.2 and abs((3.9-0.5)*(3.9+2)**2 - 118.4) < 0.2
# Kp の妥当性（1000 K）
import math
print('Kp(1000K) bar^2 ~', math.exp(-(206000-1000*215)/(8.314*1000)))

print('ALL OK')
