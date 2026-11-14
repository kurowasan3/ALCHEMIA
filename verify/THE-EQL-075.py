ITEM = '/home/claude/chem4/outh/items/THE-EQL-075/item.json'
# 検証スクリプト（グループ g07）  python3 -I この.py
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

# 小問数と答えの数の一致
n_items = it['problem'].count('\\item')
n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
print('小問数', n_items, '答え', n_ans)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from sympy import Rational as Q
RT = Q(831, 100) * 1000 * 298
H = [5, 3, 1, 0]                       # 水和数
P = [Q(1000), Q(750), Q(100)]          # 隣り合う組の解離圧（5/3, 3/1, 1/0）

def state(n_salt, n_water_total, V):
    """全水（固体＋気体）の量と体積から平衡状態を探す。戻り値 (p, {水和数: mol})"""
    # (a) 2 固体共存 (H[i], H[i+1]) で p=P[i]
    for i in range(3):
        ng = P[i] * V / RT
        ns = n_water_total - ng          # 固体中の水
        # x mol が H[i] 型、残りが H[i+1] 型
        x = (ns - n_salt * H[i+1]) / (H[i] - H[i+1])
        if 0 <= x <= n_salt:
            return P[i], {H[i]: x, H[i+1]: n_salt - x}
    # (b) 1 固体のみ
    for h in H:
        ng = n_water_total - n_salt * h
        if ng < 0: continue
        p = ng * RT / V
        # 安定範囲：h=5: p>=P0, h=3: P1<=p<=P0, h=1: P2<=p<=P1, h=0: p<=P2
        lo = {5: P[0], 3: P[1], 1: P[2], 0: 0}[h]; hi = {5: 10**9, 3: P[0], 1: P[1], 0: P[2]}[h]
        if lo <= p <= hi:
            return p, {h: n_salt}
    raise ValueError

n0 = Q(1, 100)
for V in (Q(30), Q(100)):
    p, s = state(n0, 5*n0, V)
    print('V', V, 'p', float(p), {k: sig(float(v), 2) for k, v in s.items()})
p, s = state(n0, 5*n0, Q(30)); assert p == 1000 and sig(float(s[5]), 2) == 0.0039 and sig(float(s[3]), 2) == 0.0061
p, s = state(n0, 5*n0, Q(100)); assert p == 750 and sig(float(s[3]), 2) == 0.0049 and sig(float(s[1]), 2) == 0.0051
V1 = 2*n0*RT/P[0]; V5 = 5*n0*RT/P[2]
print('V1', float(V1), 'V5', float(V5), [float(x) for x in (2*n0*RT/P[1], 4*n0*RT/P[1], 4*n0*RT/P[2])])
assert sig(float(V1), 2) == 50 and sig(float(V5), 2) == 1200
# 乾燥剤
nw = Q(2000) * 10 / RT
for n_salt, exp in ((Q(1, 100), (100, {1: 0.0077, 0: 0.0023})), (Q(4, 1000), (750, {1: 0.0035, 3: 0.00052}))):
    p, s = state(n_salt, nw, Q(10))
    got = (int(p), {k: sig(float(v), 2) for k, v in s.items()})
    print('乾燥剤', float(n_salt), got)
    assert got == exp, got
# 相対湿度
rh = Q(1000) / 3200 * 100; print('RH', float(rh)); assert round(float(rh)) == 31
assert 1000 < 0.5 * 3200 < 3100
print('すべて OK')
