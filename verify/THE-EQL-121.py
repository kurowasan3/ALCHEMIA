ITEM = '/home/claude/chemx/kyodai/outh/items/THE-EQL-121/item.json'
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
xr = float(sp.nsolve((1-x**3)/(1+x**3)-0.6, x, 0.6)); chk('xR', xr, 0.630, 3)
xs = xr**4; c = 1-(xr+xs)/2; chk('c', c, 0.61, 2)
chk('ester ee', ((1-xs)-(1-xr))/((1-xs)+(1-xr))*100, 39, 2)
c2 = 1-(0.630+0.630**4)/2; chk('c(近似値)', c2, 0.61, 2)
xr3 = (0.02/1.98)**(1/3); chk('xR98', xr3, 0.216, 3)
chk('回収', (xr3+xr3**4)/2*100, 11, 2)
chk('回収(近似)', (0.216+0.216**4)/2*100, 11, 2)
s = math.log10(0.1)/math.log10(0.9); chk('s', s, 22, 2)
chk('s(近似)', 1/(1-2*0.477), 22, 2)

print('ALL OK')
