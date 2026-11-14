ITEM = '/home/claude/chemk/outh/items/THE-ION-117/item.json'
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

from math import sqrt, log10
Ksp, KC, K1, K2, Kw, CT = 1.8e-11, 5.0e-9, 4.5e-7, 4.7e-11, 1.0e-14, 1.0e-3
OH = sqrt(Ksp/1.0e-4); pH = 14 + log10(OH); H = Kw/OH
print(OH, pH, H); assert sig(OH,2) == 4.2e-4 and round(pH,1) == 10.6
assert H/K1 < 1e-4
a2 = K2/(H+K2); a2full = K2/(H+K2+H*H/K1); print(a2, a2full); assert sig(a2,2) == 0.67 and sig(a2full,2) == 0.67
Ca = KC/(a2full*CT); print(Ca); assert sig(Ca,2) == 7.5e-6
for p, lim in [(9.0, 0.18), (10.0, 1.8e-3)]:
    oh = 10**(p-14); m = Ksp/oh**2; print(p, m); assert sig(m,2) == lim and m > 1.3e-3
H9 = 1e-9; a29 = K2/(H9+K2+H9*H9/K1); Ca9 = KC/(a29*CT); print(a29, Ca9, Ca9/Ca)
assert sig(Ca9,2) == 1.1e-4 and 14 < Ca9/Ca < 16
assert 1e-3*OH**2 < 5.5e-6   # Ca(OH)2 は沈殿しない（[Ca2+] は高々 1e-3 程度）

print('全検証OK:', ITEM)
