# 検証スクリプト（東大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/todai/out/items/THE-ELC-058/item.json'
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
F = 9.65e4
W = {'H2': 237e3, 'CO': 257e3, 'CH3OH': 702e3}
n = {'H2': 2, 'CO': 2, 'CH3OH': 6}
E = {k: W[k] / (n[k] * F) for k in W}
print(E)
assert [sig(E[k], 3) for k in ('H2', 'CO', 'CH3OH')] == [1.23, 1.33, 1.21]
assert max(E, key=E.get) == 'CO'
# メタノールの電子数：2CH3OH + 3O2 -> 2CO2 + 4H2O で O2 3 mol が 12 mol の電子を受け取る
assert check_equation('2CH3OH + 3O2 -> 2CO2 + 4H2O') is None and 12 / 2 == 6
assert sig(3.6e6 / 237e3, 2) == 15
for s in ['1.23', '1.33', '1.21', '15 mol']:
    assert s in it['answer'], s

print('ALL OK')
