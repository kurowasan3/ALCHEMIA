import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def check_all_equations(it):
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return eqs

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def check_parts(it):
    np_ = it['problem'].count('\\item')
    na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), (np_, na)
    print('小問数OK:', np_)
import sympy as sp
it = json.load(open('/home/claude/chem3/outh/items/THE-SOL-051/item.json', encoding='utf-8'))
check_all_equations(it); check_parts(it)
L = 6.0e3/18; C = 420
x1 = C*5.0/L; print('pure', x1); assert sig(x1, 2) == 6.3
print('full freeze needs', 100*L/C)
dTf = 1.86*0.0500/0.100; assert abs(dTf - 0.93) < 1e-12
x = sp.symbols('x', positive=True)
res = {}
for sc in (5.0, 8.0):
    t = -1.86*0.0500*1000/(100 - x)
    s = sp.nsolve(L*x - C*(t + sc), x, 5); tv = float(t.subs(x, s)); res[sc] = (float(s), tv)
    print(sc, float(s), tv)
    # iteration as written
    tt = -0.93
    for _ in range(4):
        xx = C*(tt + sc)/L; tt = -93/(100 - xx); print('  iter', round(xx, 3), round(tt, 4))
assert sig(res[5.0][0], 2) == 5.1 and round(res[5.0][1], 2) == -0.98
assert sig(res[8.0][0], 2) == 8.8 and round(res[8.0][1], 2) == -1.02
print('overestimate', (-res[8.0][1] - 0.93)/0.93)
print('ALL OK')
