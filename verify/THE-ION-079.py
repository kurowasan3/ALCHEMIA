# 検証スクリプト（東大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/todai/out/items/THE-ION-079/item.json'
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
from math import log10, sqrt
from sympy import symbols, nsolve
c = 2.5e-7*4.0e-7*40; s = 1.25e-5*1.2e-2*6.4e-4
assert abs(c - 4.0e-12) < 1e-20 and abs(s - 9.6e-11) < 1e-20
ph1 = -log10(sqrt(c)); ph2 = -log10(sqrt(c + s)); print(ph1, ph2)
assert round(ph1, 1) == 5.7 and round(ph2, 1) == 5.0
assert sig(s/(c+s)*100, 2) == 96
assert sig(4.0e-12/(1.25e-5*1.2e-2), 2) == 2.7e-5
assert sig(1.25e-5*1.2e-2/(2.5e-7*4.0e-7), 2) == 1.5e6
# 近似の確認：水の電離と2段目（Ka2=6e-8, CO3 は 4.7e-11）を入れる
h = symbols('h', positive=True)
Kw, Ka2s, Ka2c = 1e-14, 6e-8, 4.7e-11
eq = h - (c/h)*(1 + 2*Ka2c/h) - (s/h)*(1 + 2*Ka2s/h) - Kw/h
H = nsolve(eq, h, 1e-5); print('exact pH', -log10(float(H))); assert abs(-log10(float(H)) - ph2) < 0.01

print('ALL OK')
