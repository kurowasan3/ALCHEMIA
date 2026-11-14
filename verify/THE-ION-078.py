# 検証スクリプト（東大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/todai/out/items/THE-ION-078/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

# 小問の数と答えの数の一致
nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans == 0), (nprob, nans)
print('小問数', nprob, '答え', nans)

# ---- 個別の検証 ----
from math import sqrt, log10
K, KI = 1.8e-10, 8.5e-17
def approx(c, v, V0=50.0):
    n = c * (V0 - v) / 1000; Vt = (V0 + v) / 1000
    if abs(n) < 1e-15: return -log10(sqrt(K))
    x = abs(n) / Vt
    return -log10(x) if n > 0 else -log10(K / x)
def exact(c, v, V0=50.0):
    d_ = c * (V0 - v) / 1000 / ((V0 + v) / 1000)
    return -log10((d_ + sqrt(d_*d_ + 4*K)) / 2)
a = [approx(0.10, v) for v in (49.5, 50.0, 50.5)]; e = [exact(0.10, v) for v in (49.5, 50.0, 50.5)]
print(a, e); assert [round(x, 1) for x in a] == [3.3, 4.9, 6.4] == [round(x, 1) for x in e]
b = [approx(0.010, v) for v in (49.5, 50.5)]; eb = [exact(0.010, v) for v in (49.5, 50.5)]
print(b, eb); assert [round(x, 1) for x in b] == [4.3, 5.4]
assert (b[1] - b[0]) < (a[2] - a[0]) and (eb[1] - eb[0]) < (e[2] - e[0])
cl = 0.050 * 50 / 75; ag = K / cl; I = KI / ag; print(I); assert sig(I, 2) == 1.6e-8
assert KI / 0.050 < K / 0.050

print('ALL OK')
