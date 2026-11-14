# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-016/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and not re.search(r'_n|\bn[A-Z]|R', b):
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

data = {'銅': ('金属', True, True, None), '黒鉛': ('共有(層)', True, False, None), 'KCl': ('イオン', False, True, True),
        'HCl': ('分子', False, False, True), 'スクロース': ('分子', False, False, False), 'SiO2': ('共有', False, False, None)}
# (結晶, 固体で導電, 融解で導電, 水溶液で導電[None=水に不溶/対象外])
q1 = [k for k, v in data.items() if v[1]]
q2 = [k for k, v in data.items() if not v[1] and v[2]]
q3 = [k for k, v in data.items() if not v[1] and not v[2] and v[3] is True]
q4 = [k for k, v in data.items() if v[3] is False]
print(q1, q2, q3, q4)
assert q1 == ['銅', '黒鉛'] and q2 == ['KCl'] and q3 == ['HCl'] and q4 == ['スクロース']

print('全検証OK:', ITEM)
