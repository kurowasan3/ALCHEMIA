import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
import sympy as sp

def ce_blocks(s):
    out = []; i = 0
    while True:
        j = s.find('\\ce{', i)
        if j < 0: return out
        k = j + 4; d = 1
        while d:
            if s[k] == '{': d += 1
            elif s[k] == '}': d -= 1
            k += 1
        out.append(s[j+4:k-1]); i = k

def check_item_equations(fid):
    it = json.load(open(f'/home/claude/chem/out2/items/{fid}/item.json', encoding='utf-8'))
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for e in ce_blocks(it[key]):
            arrows = len(re.findall(r'->|<=>|<-', e))
            if arrows != 1 or '\\boxed' in e: continue
            lhs, rhs = re.split(r'<=>|->', e)
            if not lhs.strip() or not rhs.strip(): continue
            r = check_equation(e)
            assert r is None, (key, e, r)
            n += 1
    print(f'{fid}: {n} equations balanced')
    return it

def r2(x, n):  # 有効数字 n 桁の文字列
    return f'{float(x):.{n}g}'

it = check_item_equations('THE-ELC-027')
F = 96500; ne = sp.Rational(2)*4825/F
assert abs(ne - 0.1) < 1e-12
c = ne/2  # ClO- mol/L in 1.00 L
Kw, Ka = 1e-14, 3.0e-8; Kb = Kw/Ka
x = sp.symbols('x', positive=True)
exact = [s for s in sp.solve(x**2/(c - x) - Kb, x) if s > 0][0]
approx = sp.sqrt(c*Kb)
pH = 14 + sp.log(exact, 10)
print('OH exact', float(exact), 'approx', float(approx), 'ratio', float(approx/c), 'pH', float(pH))
assert round(float(pH), 1) == 10.1
assert float(approx/c) < 0.005
# with given logs
assert abs((-7 - (0.30+0.48))/2 + 14 - 10.11) < 1e-9
# (4) divided cell
pH4 = 14 + sp.log(ne/0.5, 10); print('pH4', float(pH4)); assert round(float(pH4), 1) == 13.3
# (5)
nHClO = c/(1 + Ka/1e-7); print('nHClO', float(nHClO)); assert r2(nHClO, 2) == '0.038'
print('THE-ELC-027 OK')
