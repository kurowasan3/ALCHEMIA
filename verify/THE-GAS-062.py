ITEM = '/home/claude/chemx/tokodai/outh/items/THE-GAS-062/item.json'
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
chk('(1)', float((R(123,10)+R(424,10))/2), 27.35, 4)
assert (12.3+42.4)/2 > 23.4
chk('(3)', (46.8-12.3)/42.4*100, 81.4, 3)
n = sp.symbols('n')
x0 = (R(123,10)+R(424,10))/2/1013
T = 293 + 44000*n/R(290,10)
s = sp.nsolve((x0-n)/(1-n)*1013 - (R(234,10)+R(166,100)*(T-293)), n, 0.001)
Tf = float(293+44000*s/29)
assert 293 < Tf < 298
assert abs((31.7-23.4)/5 - 1.66) < 1e-12
chk('(4)', float(s), 1.13e-3, 3)
# 1-n≈1 の近似でも同じ
s_ap = sp.solve((x0-n)*1013 - (R(234,10)+R(166,100)*(293+44000*n/29-293)), n)[0]
chk('(4)近似', float(s_ap), 1.12e-3, 3)
s2 = sp.solve((x0-n)/(1-n)*1013 - R(234,10), n)[0]
chk("(5) n'", float(s2), 3.99e-3, 3)
chk('(5)', float(s2/s), 3.5, 2)

print('ALL OK')
