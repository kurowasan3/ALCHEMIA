ITEM = '/home/claude/chemx/kyodai/out/items/THE-ION-088/item.json'
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
Ka = 2.0e-5; Kw = 1e-14; V = 0.100
h = sp.symbols('h', positive=True)
def pH(nHA, nA, nHCl):
    CT = (nHA+nA)/V; Na = nA/V; Cl = nHCl/V
    f = lambda p: (10**-p) + Na - CT*Ka/(Ka+10**-p) - Cl - Kw/10**-p
    lo, hi = 0.0, 14.0
    for _ in range(200):
        mid = (lo+hi)/2
        if f(mid) > 0: lo = mid
        else: hi = mid
    return (lo+hi)/2
chk('P', pH(1e-2,1e-2,0), 4.7, 2); chk('Q', pH(1e-3,1e-3,0), 4.7, 2)
chk('P+HCl', pH(1e-2,1e-2,2e-3), 4.5, 2); chk('Q+HCl', pH(1e-3,1e-3,2e-3), 2.0, 2)
chk('Pの近似', 4.699+0.301-0.477, 4.5, 2)
chk('倍率', (9/11*1e-2)/(9/11*1e-3), 10, 2)
dP = pH(1e-2,1e-2,0)-pH(1e-2,1e-2,4e-3); dR = pH(4e-3,1.6e-2,0)-pH(4e-3,1.6e-2,4e-3)
print('厳密な変化', dP, dR); assert dP < dR
chk('dP', dP, 0.37, 2); chk('dR', dR, 0.43, 2); chk('dP近似', 0.845-0.477, 0.37, 2); chk('dR近似', 0.602-(0.477-0.301), 0.43, 2)

print('ALL OK')
