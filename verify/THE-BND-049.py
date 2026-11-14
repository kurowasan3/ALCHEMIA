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
import sympy as sp, math
it = json.load(open('/home/claude/chem3/outh/items/THE-BND-049/item.json', encoding='utf-8'))
check_all_equations(it); check_parts(it)
A1, A2 = 1.748, 1.763
rho = sp.symbols('rho', positive=True)
for s3 in (sp.sqrt(3), 1.732):
    sol = sp.solve(sp.Eq(A2/s3, A1/(1 + rho)), rho); print('crossover', [float(v) for v in sol])
    assert round(float(sol[0]), 3) == 0.717
c1, c2 = math.sqrt(2) - 1, math.sqrt(3) - 1
print('critical', c1, c2); assert round(c1, 3) == 0.414 and round(c2, 3) == 0.732
def f(r):
    dN = 1 + r if r >= c1 else math.sqrt(2); dC = 1 + r if r >= c2 else math.sqrt(3)
    return A1/dN, A2/dC
import numpy as np
for r in np.linspace(0.01, 1.0, 991):
    n, c = f(r)
    assert (c > n) == (r > 0.71731), r
print('rho=0.72', f(0.72)); print('diff %', (A2 - A1)/A1*100); assert sig((A2 - A1)/A1*100, 2) == 0.86
assert round(1.732/1.72, 3) == 1.007
print('ALL OK')
