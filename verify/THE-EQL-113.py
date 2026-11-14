ITEM = '/home/claude/chemx/kyodai/out/items/THE-EQL-113/item.json'
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
x = sp.symbols('x', positive=True)
K = 1e-3; C = 2e-3
a = [s for s in sp.solve(x**2/(C-x)-K, x) if s < C][0]; chk('α', float(a)/C, 0.50, 2)
z = [s for s in sp.solve((x/5)*x/(C-x)-K, x) if s < C][0]; chk('β', float(z)/C, 0.77, 2)
chk('β(√65=8.06)', (-5.000+8.062)/2/2, 0.77, 2)
chk('吸着割合', 0.4/0.5, 0.80, 2)
A2, P2, R2 = C-float(z), float(z)/5, float(z)
assert A2 < float(a) and P2 < float(a) and R2 > float(a)
w = [s for s in sp.solve(x**2/((C-x)/5)-K, x) if s < C][0]; chk('A吸着の反応率', float(w)/C, 0.27, 2)
assert float(w)/C < 0.5

print('ALL OK')
