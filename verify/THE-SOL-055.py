# 検証スクリプト（グループ g03）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, N
ITEM = '/home/claude/chem4/out/items/THE-SOL-055/item.json'
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

for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1:
            r = check_equation(b); assert r is None, (b, r); print('式OK:', b)

def sig(x, n=3):
    from math import log10, floor
    x = float(x); return round(x, -int(floor(log10(abs(x)))) + (n - 1))

mM = Q('40.0') * Q('0.791')          # メタノールの質量
mW = Q('63.3') - mM                  # 水の質量
n = mM / Q('32.0')
assert sig(mM / Q('63.3') * 100) == 50.0
assert sig(n / (mW / 1000)) == 31.2
V = Q('63.3') / Q('0.916')           # mL
assert sig(V) == 69.1
assert sig(n / (V / 1000)) == 14.3
assert sig(Q('40.0') / V * 100) == 57.9
Vadd = Q('40.0') + mW / Q('0.998')
assert sig(Vadd) == 71.7
assert sig(n / (Vadd / 1000)) == 13.8 and n / (Vadd / 1000) < n / (V / 1000)
print('収縮量 mL =', N(Vadd - V, 3))
for s in ('50.0\\,\\%', '31.2\\,\\mathrm{mol/kg}', '14.3\\,\\mathrm{mol/L}', '小さくなる'):
    assert s in it['answer'], s
print('THE-SOL-055 OK')
