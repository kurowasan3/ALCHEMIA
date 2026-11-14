ITEM = '/home/claude/chemx/kyodai/outh/items/THE-EQL-114/item.json'
# 検証スクリプト（京大 g01）  python3 -I この.py
import sys, json, re, math
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
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))

def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)) or v == expect, (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')

# 小問の数と答えの数の一致
nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), ('小問数', nq, na)
print('小問数OK', nq)

# ---- 個別の検証 ----

import sympy as sp
t, k, M, rho, r0 = sp.symbols('t k M rho r0', positive=True)
r = sp.Function('r')
sol = sp.dsolve(sp.Eq(sp.diff(rho/M*sp.Rational(4,3)*sp.pi*r(t)**3, t), -k*4*sp.pi*r(t)**2), r(t), ics={r(0): r0})
sols = sol if isinstance(sol, list) else [sol]
rt = [s.rhs for s in sols if sp.simplify(s.rhs.subs(t, 0) - r0) == 0][0]
print('r(t) =', sp.simplify(rt))
assert sp.simplify(sp.diff(rt, t, 2)) == 0
tau = sp.solve(sp.Eq(rt, 0), t)[0]; assert sp.simplify(tau - rho*r0/(k*M)) == 0
n0 = rho/M*sp.Rational(4,3)*sp.pi*r0**3; v0 = k*4*sp.pi*r0**2
assert sp.simplify(tau/(n0/v0)) == 3
chk('半分の時刻', 1-0.5**(1/3), 0.206, 3); chk('半分の時刻(近似)', 1-0.794, 0.206, 3)
chk('速さ', 0.5**(2/3), 0.630, 3); chk('速さ(近似)', 0.794**2, 0.630, 3)
chk('tau', 3*(0.280/56)/1.00e-6, 1.50e4, 3); chk('tau/10', 3*(0.280/56)/1.00e-6/10, 1.50e3, 3)
chk('混合', (8*(1-(1/2)**3)+1)/9*100, 88.9, 3)

print('ALL OK')
