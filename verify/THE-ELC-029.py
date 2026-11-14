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

it = check_item_equations('THE-ELC-029')
F = 96500
ne = sp.Rational(77200, 96500); assert ne == sp.Rational(4, 5)
H2 = sp.Rational(2, 100); Fe3 = ne; Cr2 = ne - 2*H2
assert Fe3 == sp.Rational(800, 1000) and Cr2 == sp.Rational(760, 1000)
dis = min(Fe3, Cr2); Q = dis*F; print('Q', float(Q), 'eff', float(dis/ne))
assert r2(Q, 3) == '7.33e+04' and dis/ne == sp.Rational(95, 100)
left = Fe3 - dis; assert left == sp.Rational(4, 100)
assert 2*H2 == left  # H2 exactly reduces leftover Fe3+
a, b = sp.symbols('a b')
sol = sp.solve([Cr2 - a - b - sp.Rational(700, 1000), a - sp.Rational(3, 100)], [a, b]); print(sol)
assert sol[b] == sp.Rational(3, 100)
Fe3_after = Fe3 - sol[a] - sol[b] - sp.Rational(700, 1000); assert Fe3_after == sp.Rational(4, 100)
print('THE-ELC-029 OK')
