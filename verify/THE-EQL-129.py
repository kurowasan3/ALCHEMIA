ITEM = '/home/claude/chemx/tokodai/outh/items/THE-EQL-129/item.json'
# 検証スクリプト（東工大 g04）  python3 -I この.py
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
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')

nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na, ('小問数', nq, na)
print('小問数OK', nq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
# ---- 個別の検証 ----

import sympy as sp
R = sp.Rational
K = R(15,10)*R(15,10)/(R(15,10)*R(5,10))
assert K == 3
y = sp.symbols('y')
ys = sp.solve(283*y + 242*(R(8,10)-y) - R(2059,10), y)[0]
chk('(2)CO', float(ys), 0.300, 3); chk('(2)H2', float(R(8,10)-ys), 0.500, 3)
def final(yco):
    zh = R(8,10)-yco
    co, h2o, co2, h2 = R(15,10)-yco, R(5,10)+zh, R(15,10)+yco, R(15,10)-zh
    x = sp.symbols('x')
    sols = [s for s in sp.solve((co2+x)*(h2+x) - 3*(co-x)*(h2o-x), x) if all(v > 0 for v in (co-s, h2o-s, co2+s, h2+s))]
    assert len(sols) == 1
    s = sols[0]
    return (co2*h2/(co*h2o), s, co-s, h2o-s, co2+s, h2+s)
Q, x, co, h2o, co2, h2 = final(ys)
chk('Q', float(Q), 1.5, 2); chk('x', float(x), 0.200, 3)
xx=sp.symbols('xx'); assert set(sp.solve(2*xx**2-R(94,10)*xx+R(18,10),xx))=={R(1,5),R(9,2)}
assert sp.expand((R(18,10)+xx)*(1+xx)-3*(R(12,10)-xx)*(1-xx)) == sp.expand(-(2*xx**2-R(94,10)*xx+R(18,10)))
for nm,v,e in [('CO',co,1.00),('H2O',h2o,0.800),('CO2',co2,2.00),('H2',h2,1.20)]:
    chk(nm, float(v), e, 3)
for alt in (0, R(8,10)):
    r = final(alt)
    assert sp.simplify(r[5]-h2) == 0, alt
h = sp.symbols('h')
assert R(6,5) in sp.solve(h*(R(4,5)+h) - 3*(2-h)*(R(11,5)-h), h)
chk('(5)', float(283*(R(15,10)-co) + 242*(R(15,10)-h2)), 214.1, 4)
assert R(2059,10) + 41*x == 283*(R(15,10)-co) + 242*(R(15,10)-h2)
for e in ['CO + H2O <=> CO2 + H2','2CO + O2 -> 2CO2','2H2 + O2 -> 2H2O']:
    assert check_equation(e) is None

print('ALL OK')
