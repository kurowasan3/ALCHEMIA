ITEM = '/home/claude/chemk/out/items/THE-ELC-088/item.json'
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

k1 = 120*0.10/1000; R = 5.0/(k1*4.0); I = 4.2/R; Q = I*1800; ne = Q/9.65e4; gas = ne*3/4; V = gas*22.4e3
print(k1, R, I, Q, ne, gas, V)
assert sig(R,3) == 104 and sig(I,3) == 4.03e-2 and sig(V,3) == 12.6
k2 = 113*0.20/1000; r = k2/k1; print(r); assert sig(r,3) == 1.88
assert check_equation('2H2O -> 2H2 + O2') is None

print('全検証OK:', ITEM)
