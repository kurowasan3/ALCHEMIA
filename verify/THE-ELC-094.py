ITEM = '/home/claude/chemk/outh/items/THE-ELC-094/item.json'
# 検証スクリプト（グループ g04）  python3 -I この.py
import sys, json, re
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
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from math import log10
S = 134.0-76.0; E0 = 76.0 - S*(-3); print(S, E0); assert S == 58.0 and E0 == 250.0
k = 10**((35.5-E0)/S)/0.100; print('k', k); assert sig(k,2) == 2.0e-3
Ef = E0 + S*log10(4.0e-3 + k*0.140); print('順方向 E', Ef); assert abs(Ef-112.6) < 0.05
naive = 10**((112.6-E0)/S); corr = naive - 2.0e-3*0.140
print(naive, corr, (naive-corr)/corr*100)
assert sig(naive,2) == 4.3e-3 and sig(corr,2) == 4.0e-3 and sig((naive-corr)/corr*100,2) == 7.0
kmax = 0.010*4.0e-3/0.140; print(kmax); assert sig(kmax,2) == 2.9e-4
assert sig(2.0e-3*0.140,2) == 2.8e-4
# 希釈：比が不変
E10 = E0 + S*log10((4.0e-3 + k*0.140)/10); assert abs((Ef-E10)-S) < 1e-9
naive10 = 10**((E10-E0)/S); corr10 = naive10 - k*0.014
assert abs((naive10-corr10)/corr10 - (k*0.140)/4.0e-3) < 1e-12

print('全検証OK:', ITEM)
