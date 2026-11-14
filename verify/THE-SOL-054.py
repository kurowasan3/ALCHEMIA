# 検証スクリプト（グループ g03）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q
ITEM = '/home/claude/chem4/out/items/THE-SOL-054/item.json'
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

neq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1:
            r = check_equation(b); assert r is None, (b, r); print('式OK:', b); neq += 1
assert neq == 2

def sig(x, n=2):
    from math import log10, floor
    x = float(x); return round(x, -int(floor(log10(abs(x)))) + (n - 1))

glc = Q('0.10'); cacl2 = Q('0.040') * 3
urea = Q(1000) * Q('1.00') * Q('0.60') / 100 / (12 + 16 + 2 * (14 + 2 * Q('1.0')))
nacl = Q('0.10') * 2; na2so4 = Q('0.070') * 3
assert cacl2 > glc                    # (1) 右が高い
assert urea == glc                    # (2) 変化なし
assert na2so4 > nacl                  # (3) 右が高い
R = Q(83, 10) * 1000
d = (cacl2 - glc) * R * 300
assert sig(d) == 5.0e4                # (4)
assert (cacl2 - glc) * R * 320 / d == Q(320, 300)  # (5) T に比例
print('THE-SOL-054 OK')
