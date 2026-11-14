ITEM = '/home/claude/chem4/outh/items/THE-ION-067/item.json'
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

from sympy import Rational as Q, symbols, solve
K = Q(10, 10**11) / Q(50, 10**10)
print('K', K); assert K == Q(1, 50)

def treat(n_solid_A, n_anion, K):
    """固体 A（変換される側）n_solid_A と陰イオン n_anion（体積共通）。A が残ると仮定して x を求め、超えれば全量変換"""
    x = K * n_anion / (1 + K)
    if x < n_solid_A:
        return x, True          # A が残る
    return n_solid_A, False     # 全量変換
nA = Q(5, 1000); nC = Q(1, 10)
x, remain = treat(nA, nC, K)
print('1回目', float(x), remain); assert remain and sig(float(x), 3) == 1.96e-3
Ba = Q(50, 10**10) / ((nC - x) / Q(1, 10))
print('[Ba2+]', float(Ba)); assert sig(float(Ba), 2) == 5.1e-9
# くり返し
left = nA; k = 0
while left > 0:
    k += 1
    x, remain = treat(left, nC, K)
    left -= x
    print(k, '回目 変換', float(x), '残り', float(left))
assert k == 3
SO4 = x / Q(1, 10); CO3 = (nC - x) / Q(1, 10)
print('[SO4]', float(SO4), '比', float(SO4/CO3)); assert sig(float(SO4), 3) == 1.08e-2 and SO4/CO3 < K
# 逆向き
y, remain = treat(Q(5, 1000), Q(1, 100), 1/K)
print('逆向き', float(y), remain); assert not remain and y == Q(5, 1000)
yy = (1/K) * Q(1, 100) / (1 + 1/K); print('仮定時の y', float(yy)); assert yy > Q(5, 1000)
c_CO3 = Q(5, 1000) / Q(1, 10); c_SO4 = (Q(1, 100) - Q(5, 1000)) / Q(1, 10)
assert sig(float(c_CO3), 2) == 0.05 and sig(float(c_SO4), 2) == 0.05 and c_CO3 / c_SO4 < 1/K
# (5) 2.00 mol/L × 50.0 mL
x2, _ = treat(nA, Q(2) * Q(50, 1000), K)
assert x2 == treat(nA, nC, K)[0]
print('すべて OK')
