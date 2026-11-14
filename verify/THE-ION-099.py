# 検証スクリプト THE-ION-099  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/outh/items/THE-ION-099/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    if n_items: assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()


import numpy as np
def comp(pH, p1, p2):
    r = np.array([1, 10**(pH-p1), 10**(2*pH-p1-p2)]); r = r/r.sum()
    return r, -(r[1] + 2*r[2])
PH, TP = (2.9, 5.4), (3.5, 4.8)
for p in (PH, TP): assert abs(comp(4.15, *p)[1] + 1) < 1e-12
vals = {pH: (round(comp(pH,*PH)[1],2), round(comp(pH,*TP)[1],2)) for pH in (2.0,3.2,4.15,5.1,7.0)}
print(vals)
assert vals[3.2] == (-0.67, -0.35) and vals[5.1] == (-1.33, -1.65)
assert vals[2.0] == (-0.11, -0.03) and vals[7.0] == (-1.98, -1.99)
d = {pH: abs(comp(pH,*PH)[1]-comp(pH,*TP)[1]) for pH in vals}
best = max(d.values()); print(d)
assert abs(d[3.2]-d[5.1]) < 1e-9 and abs(d[3.2]-best) < 1e-12
x = np.linspace(1, 9, 80001); dz = np.array([comp(v,*PH)[1]-comp(v,*TP)[1] for v in x])
print('極値の pH', x[dz.argmin()], x[dz.argmax()])
# 近似値での手計算
z = lambda a, b, c: -(b + 2*c)/(a + b + c)
assert round(z(1,2.0,0.0126),2) == -0.67 and round(z(1,0.50,0.0126),2) == -0.35
assert round(z(0.0063,1,0.50),2) == -1.33 and round(z(0.025,1,2.0),2) == -1.65
# 中点での HA- の比率と傾き
print('HA-', comp(4.15,*PH)[0][1], comp(4.15,*TP)[0][1])
e = 1e-5
s = [(comp(4.15+e,*p)[1]-comp(4.15-e,*p)[1])/(2*e) for p in (PH, TP)]
print('傾き', s); assert abs(s[1]) > abs(s[0])
print('OK')
