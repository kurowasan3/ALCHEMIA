ITEM = '/home/claude/chemx/tokodai/out/items/THE-GAS-060/item.json'
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
p0 = sp.Rational(6,10)*sp.Rational(234,10)
chk('(1)', float(p0*303/293/sp.Rational(424,10)*100), 34.2, 3)
chk('(2)', float(p0/sp.Rational(424,10)*100), 33.1, 3)
T = sp.symbols('T')
sol = sp.solve(sp.Eq(p0*T/293, sp.Rational(123,10)+sp.Rational(94,100)*(T-283)), T)[0]
assert 283 < sol < 288
assert float(p0*283/293) > 12.3 and float(p0*288/293) < 17.0
chk('(3)', float(sol), 284.4, 4)
chk('(4)', float(sp.Rational(123,10)*293/283/sp.Rational(234,10)*100), 54.4, 3)
# 区間の傾き
assert abs((17.0-12.3)/5 - 0.94) < 1e-12

print('ALL OK')
