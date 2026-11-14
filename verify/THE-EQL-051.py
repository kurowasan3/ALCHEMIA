# 検証スクリプト THE-EQL-051  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-EQL-051/item.json'
import json, re
from sympy import symbols, solve, simplify, limit, oo, Rational as R
it = json.load(open(ITEM, encoding='utf-8'))
n_items = len(re.findall(r'\\item\b', it['problem'])); n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans == 4, (n_items, n_ans)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
k1, k2, k3, E0, S, ES = symbols('k1 k2 k3 E0 S ES', positive=True)
sol = solve(k1*(E0 - ES)*S - (k2 + k3)*ES, ES)[0]
K = (k2 + k3)/k1
assert simplify(sol - E0*S/(K + S)) == 0
v = k3*sol
assert simplify(limit(v/S, S, 0) - k3*E0/K) == 0   # 1 次
assert simplify(limit(v, S, oo) - k3*E0) == 0      # 0 次
f = lambda s, Kv=R(4, 10**4), e0=R(2, 10**8), kk=50: kk*e0*s/(Kv + s)
v1 = f(R(4, 10**4)); print('v', float(v1)); assert sig(float(v1), 2) == 5.0e-7
x = symbols('x', positive=True); s90 = solve(x/(R(4, 10**4) + x) - R(9, 10), x)[0]; print('S90', s90); assert sig(float(s90), 2) == 3.6e-3
ratio_a = f(R(8, 100)) / f(R(4, 100)); print('(a)', float(ratio_a)); assert abs(float(ratio_a) - 1.005) < 0.001
ratio_b = f(R(4, 100), e0=R(4, 10**8)) / f(R(4, 100)); assert ratio_b == 2
assert abs(float(R(100, 101)) - 0.990) < 0.001 and abs(float(R(200, 201)) - 0.995) < 0.001
print('全検証OK:', ITEM)
# 反応式：validate と同じ規則（boxed・\text・\mathrm・$ を含む式は対象外）で，検査対象の式が無いことを確認
def ce_blocks(s):
    return re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', s)
eqs = [b for k in ('problem', 'answer', 'solution') for b in ce_blocks(it[k])
       if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm|\$', b)]
print('検査対象の式', eqs); assert eqs == []
