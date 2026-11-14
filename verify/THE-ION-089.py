ITEM = '/home/claude/chemx/kyodai/out/items/THE-ION-089/item.json'
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

K1 = 4.5e-7; K2 = 4.7e-11; Kw = 1e-14
C = 3.0e-7*40; chk('C', C, 1.2e-5, 2)
def pH(no3, exact=False):
    lo, hi = 0.0, 14.0
    for _ in range(200):
        p = (lo+hi)/2; h = 10**-p
        f = h - K1*C/h - no3 - ((2*K1*K2*C/h**2 + Kw/h) if exact else 0)
        if f > 0: lo = p
        else: hi = p
    return (lo+hi)/2
chk('pH0', pH(0), 5.6, 2); assert abs(pH(0)-pH(0,True)) < 0.01
h = 10**-4.5; no3 = h - K1*C/h; chk('NO3', no3, 3.1e-5, 2)
chk('寄与', (K1*C/h)/h*100, 0.54, 2)
p5 = pH(2.0e-6); chk('pH5', p5, 5.5, 2); assert abs(p5-pH(2.0e-6,True)) < 0.01
chk('x', (2.0e-6+(4.0e-12+21.6e-12)**0.5)/2, 3.5e-6, 2)

print('ALL OK')
