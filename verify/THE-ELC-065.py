ITEM = '/home/claude/chemx/kyodai/outh/items/THE-ELC-065/item.json'
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
x = sp.symbols('x')
E = 0.960/80; chk('電子', E, 1.20e-2, 3)
nCu = 2.00e-3
x1 = sp.solve(sp.Eq(sp.Rational(3,4)*(E-x), 1.80*x/4), x)[0]; chk('仮定1', float(x1), 7.50e-3, 3)
assert x1 > 2*nCu
x2 = sp.solve(sp.Eq(sp.Rational(3,4)*(E-x), 1.80*(x/4 + (x-2*nCu)/2)), x)[0]; chk('x', float(x2), 6.00e-3, 3)
assert x2 > 2*nCu
chk('y', float(E-x2), 6.00e-3, 3)
g1 = float(x2/4 + (x2-2*nCu)/2); g2 = 0.75*float(E-x2); chk('比', g2/g1, 1.80, 3)
chk('H2 mL', float((x2-2*nCu)/2)*22.4e3, 22.4, 3)
chk('正極', E*32, 0.384, 3)
assert check_equation('Pb + PbO2 + 2H2SO4 -> 2PbSO4 + 2H2O') is None

print('ALL OK')
