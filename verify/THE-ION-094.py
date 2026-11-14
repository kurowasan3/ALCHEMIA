ITEM = '/home/claude/chemx/kyodai/outh/items/THE-ION-094/item.json'
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

Kc, Ks1, Ks2, Kc2, Kw = 4.5e-7, 1.5e-2, 6.0e-8, 4.7e-11, 1e-14
C = 3.0e-7*40; S = 1.2e-5*1.0e-3
def H(S, C, exact):
    lo, hi = 0.0, 14.0
    for _ in range(200):
        p = (lo+hi)/2; h = 10**-p
        f = h - Ks1*S/h*(1 + (2*Ks2/h if exact else 0)) - Kc*C/h*(1 + (2*Kc2/h if exact else 0)) - (Kw/h if exact else 0)
        if f > 0: lo = p
        else: hi = p
    return 10**-((lo+hi)/2)
h = H(S, C, False); he = H(S, C, True)
chk('pH', -math.log10(h), 4.9, 2); assert abs(math.log10(h/he)) < 0.01
chk('SO3/HSO3', Ks2/h, 4.4e-3, 2)
chk('比', Ks1*S/(Kc*C), 33, 2)
chk('全S', S*(1+Ks1/h), 1.3e-5, 2); chk('全S pH3', S*(1+Ks1/1e-3), 1.9e-7, 2)
p = (1e-8 - Kc*C)/(Ks1*1.2e-5); chk('p', p, 5.6e-2, 2)
h4 = H(1.2e-5*p, C, True); assert abs(-math.log10(h4) - 4.0) < 0.01

print('ALL OK')
