# 検証スクリプト（グループ g03）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q
ITEM = '/home/claude/chem4/outh/items/THE-SOL-057/item.json'
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

def sig(x, n=2):
    from math import log10, floor
    x = float(x); return round(x, -int(floor(log10(abs(x)))) + (n - 1))

# 塩: (凝析価, 陽イオン(価数,個数), 陰イオン(価数,個数))
data = {
    'NaCl':  (Q(52),        (1, 1), (1, 1)),
    'KNO3':  (Q(51),        (1, 1), (1, 1)),
    'Na2SO4': (Q(26),       (1, 2), (2, 1)),
    'MgCl2': (Q('0.80'),    (2, 1), (1, 2)),
    'MgSO4': (Q('0.81'),    (2, 1), (2, 1)),
    'AlCl3': (Q('0.072'),   (3, 1), (1, 3)),
}
def czeta(salt, sign):
    c, cat, an = data[salt]
    z, nu = cat if sign == 'neg' else an   # 負のコロイド → 陽イオンが対イオン
    return c * nu * z**6
neg = [czeta(s, 'neg') for s in data]
pos = [czeta(s, 'pos') for s in data]
spread = lambda v: float(max(v) / min(v))
print('負と仮定した c z^6:', [float(x) for x in neg], 'ばらつき', spread(neg))
print('正と仮定した c z^6:', [float(x) for x in pos], 'ばらつき', spread(pos))
assert spread(neg) < 1.05 and spread(pos) > 10       # 符号は負に一意
assert [sig(czeta(s, 'neg')) for s in ('Na2SO4', 'MgCl2', 'AlCl3')] == [52, 51, 52]
c1 = Q(52)
assert sig(c1 / 3) == 17                 # K3PO4
assert sig(c1 / 3**6 / 2) == 0.036       # Al2(SO4)3
assert sig((c1 - 30) / 2**6) == 0.34     # 混合
# Y: NaCl 9.6，Na2SO4 が約 1/60 → 正（負なら 1/2）
cY = Q('9.6'); assert abs(float((cY / 64) / cY) - 1/60) < 0.003
assert sig(cY / 2) == 4.8
for s in ('負に帯電', '17\\,\\mathrm{mmol/L}', '0.036\\,\\mathrm{mmol/L}', '0.34\\,\\mathrm{mmol/L}', '4.8\\,\\mathrm{mmol/L}', '正に帯電'):
    assert s in it['answer'], s
print('THE-SOL-057 OK')
