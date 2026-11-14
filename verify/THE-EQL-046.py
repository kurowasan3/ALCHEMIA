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

it = check_item_equations('THE-EQL-046')
a, b, z = sp.symbols('a b z', real=True)
Kp = 1; sols = []
# case (i): a <= 2b
for s in sp.solve([a + b - 2, b + a/2 - z - sp.Rational(5, 4), z - Kp*(a - 2*z)**2], [a, b, z], dict=True):
    ok = s[a] <= 2*s[b] and s[a] - 2*s[z] >= 0 and s[z] >= 0 and s[a] > 0 and s[b] > 0
    print('i', s, ok)
    if ok: sols.append((s[a], s[b], s[z], 'O2', s[b] - s[a]/2))
for s in sp.solve([a + b - 2, a - z - sp.Rational(5, 4), z - Kp*(2*b - 2*z)**2], [a, b, z], dict=True):
    ok = s[a] > 2*s[b] and 2*s[b] - 2*s[z] >= 0 and s[z] >= 0 and s[b] > 0
    print('ii', s, ok)
    if ok: sols.append((s[a], s[b], s[z], 'NO', s[a] - 2*s[b]))
print(sols)
assert sorted((x[0], x[1]) for x in sols) == [(1, 1), (sp.Rational(3, 2), sp.Rational(1, 2))]
assert all(x[2] == sp.Rational(1, 4) and x[4] == sp.Rational(1, 2) for x in sols)
# total pressure as function of a is minimal at a = 4/3
def P(av):
    av = sp.nsimplify(av); bv = 2 - av; zz = sp.symbols('zz', positive=True)
    if av <= 2*bv:
        zr = [r for r in sp.solve(zz - (av - 2*zz)**2, zz) if av - 2*r >= 0][0]; return bv + av/2 - zr
    zr = [r for r in sp.solve(zz - (2*bv - 2*zz)**2, zz) if 2*bv - 2*r >= 0][0]; return av - zr
vals = [(x, float(P(x))) for x in [0.8, 1.0, 1.2, sp.Rational(4, 3), 1.4, 1.5, 1.8]]
print(vals); m = min(vals, key=lambda t: t[1]); assert m[0] == sp.Rational(4, 3)
print('THE-EQL-046 OK')
