# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import sqrt
from sympy import symbols, solve, Eq, nsimplify
ITEM = '/home/claude/chem2/outh/items/THE-ION-048/item.json'
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
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and ' + ' in b:
            r = check_equation(b); assert r is None, (b, r); print('式OK:', b)

L0 = 425 + 101 - 126; print('(1)', L0); assert L0 == 400
data = [(4.0e-3, 80.0), (4.0e-4, 200.0)]
Kas = []
for C, L in data:
    a = L / L0; Ka = C*a*a/(1 - a); Kas.append(Ka); print('(2)', C, a, Ka)
assert all(abs(k - 2.0e-4) < 1e-12 for k in Kas)
X = symbols('X', positive=True)
sols = solve(Eq(nsimplify('4.0e-3')*80**2/(X*(X - 80)), nsimplify('4.0e-4')*200**2/(X*(X - 200))), X)
sols = [s for s in sols if s > 200]; print('(3) Λ0 =', sols); assert len(sols) == 1 and sols[0] == 400
print('(3) Ka =', 25.6/(400*320)); assert abs(25.6/(400*320) - 2.0e-4) < 1e-12
print('(4)', sqrt(10), 200/80); assert round(sqrt(10), 1) == 3.2 and 200/80 == 2.5
Ka = 2.0e-4; C = 1.0e-5
a = (-Ka + sqrt(Ka*Ka + 4*Ka*C))/(2*C); print('(5) 厳密 α', a, 'Λ', a*400, ' 与えた値で', (2.19 - 2.00)*1e-4/(2.0e-5)*400)
assert round(a*400, -1) == 380 and round((2.19 - 2.00)*1e-4/(2.0e-5)*400, -1) == 380
# Λ0 に十分近づく濃度：α>0.99 となる C
Cq = Ka*(1 - 0.99)/0.99**2; print('α=0.99 となる C', Cq); assert Cq < 1e-5
print('全検証OK:', ITEM)
