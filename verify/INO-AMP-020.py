# 検証スクリプト（東大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/todai/out/items/INO-AMP-020/item.json'
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
nAl = 2.70/27; nK = 0.150
nH2 = nAl*3/2; assert sig(nH2*22.4, 3) == 3.36
left = nK - nAl
v1 = left/2/3.00*1000; v2 = (left + nAl + 3*nAl)/2/3.00*1000; print(v1, v2)
assert sig(v1, 3) == 8.33 and sig(v2, 3) == 75.0
nSO4 = (left + 4*nAl)/2
M = 27 + 39 + 2*96 + 12*18; assert M == 474
assert min(nK, nAl, nSO4/2) == nAl and sig(nAl*M, 3) == 47.4
assert check_equation('[Al(H2O)6]^3+ <=> [Al(OH)(H2O)5]^2+ + H+') is None

print('ALL OK')
