# 検証スクリプト（グループ g03）  python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-EQL-152/item.json'
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
T, w = sp.symbols('T w')
fc = lambda T: 18.0e3/T - 32.0
fy = lambda T: 13.0e3/T - 17.0
# データの作り方：真の T=293 K, δw=-1.0
dc = -1.0 + fc(293.0); dy = -1.0 + fy(293.0)
print('δc, δY (真値から)', dc, dy)
assert round(dc, 2) == 28.43 and round(dy, 2) == 26.37
# (1)
a, b = fc(273.0), fc(303.0); print(a, b)
assert round(a, 1) == 33.9 and round(b, 1) == 27.4 and a > b
# (2)
T2 = 18.0e3/(28.43 + 32.0); print('T2', T2); assert sig(T2, 3) == 298
# (3)
sol = sp.solve([sp.Eq(28.43 - w, 18.0e3/T - 32.0), sp.Eq(26.37 - w, 13.0e3/T - 17.0)], [T, w], dict=True)[0]
print(sol); assert sig(float(sol[T]), 3) == 293 and round(float(sol[w]), 1) == -1.0
assert abs((T2 - float(sol[T])) - 4.8) < 0.3
# (4) 真の δw>0 で δw=0 と仮定すると T は低く出る
Ttrue = 285.0; wtrue = 1.2; dcg = wtrue + fc(Ttrue)
Test = 18.0e3/(dcg + 32.0); print('氷期', Ttrue, Test); assert Test < Ttrue
# (5) 1/T の係数が同じなら差は定数
assert abs((fc(280.0) - (18.0e3/280.0 - 25.0)) - (-7.0)) < 1e-9 and abs((fc(310.0) - (18.0e3/310.0 - 25.0)) - (-7.0)) < 1e-9
print('OK')
