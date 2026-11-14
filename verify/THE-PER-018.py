# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem4/out/items/THE-PER-018/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        if ' + ' not in b and not re.search(r'(^|\s)\d+\s*[A-Z(\[]', b): continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

NOBLE = [2, 10, 18, 36, 54, 86, 118]
def pos(z):
    prev = 0
    for per, g in enumerate(NOBLE, start=1):
        if z <= g: break
        prev = g
    else:
        per = 8; prev = 118
    d = z - prev
    if per in (6, 7):
        if d >= 17: d -= 14
        elif d > 3: d = 3   # ランタノイド・アクチノイド
    return per, d
known = {38: (5, 2), 52: (5, 16), 79: (6, 11), 114: (7, 14), 118: (7, 18), 87: (7, 1), 88: (7, 2), 72: (6, 4), 82: (6, 14), 26: (4, 8), 119: (8, 1)}
for z, v in known.items():
    assert pos(z) == v, (z, pos(z), v)
    print(z, '→ 第%d周期 %d族' % pos(z))
assert sum([2, 8, 18, 32, 32, 18, 4]) == 114
assert check_equation('2Fr + 2H2O -> 2FrOH + H2') is None
print('OK')
