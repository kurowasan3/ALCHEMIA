# 検証スクリプト（東工大 g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/out/items/THE-SOL-089/item.json'
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

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans <= 1), (nprob, nans)
print('小問数', nprob, '答え', nans)
for k in ['src_no', 'src_label', 'field', 'theme', 'difficulty', 'time_min', 'format', 'problem', 'answer', 'solution', 'core', 'tags', 'verified', 'verify_note', 'problem_fig', 'solution_fig']:
    assert k in it, k
# ---- 個別の検証 ----

import sympy as sp
T = sp.symbols('T')
y = sp.Rational(250, 10000)*T - sp.Rational(865, 100)
x = sp.Rational(400, 10000)*T - sp.Rational(1360, 100)
Td = sp.solve(sp.Eq(y, sp.Rational(40, 100)), T)[0]; print('dew', float(Td)); assert round(float(Td)) == 362
x1 = x.subs(T, Td); print('first drop', float(x1)); assert abs(float(x1) - 0.88) < 1e-9
L = (0.40 - float(y.subs(T, 356))) / (float(x.subs(T, 356)) - float(y.subs(T, 356))); print('L', L); assert round(L, 2) == 0.38
Tb = sp.solve(sp.Eq(x, sp.Rational(40, 100)), T)[0]; print('bubble', float(Tb)); assert float(Tb) == 350
assert abs(float(y.subs(T, Tb)) - 0.10) < 1e-9
T5 = sp.solve(sp.Eq(y, sp.Rational(30, 100)), T)[0]; print('dew z=0.30', float(T5)); assert float(T5) == 358
for t in range(347, 366):
    xv, yv = float(x.subs(T, t)), float(y.subs(T, t))
    assert 0 <= yv < xv <= 1, (t, xv, yv)
print('OK')
