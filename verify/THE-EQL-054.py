# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, Rational as Q, sqrt, solve, nsolve, diff, N
ITEM = '/home/claude/chem2/outh/items/THE-EQL-054/item.json'
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; d = 1
        while d:
            if s[j] == '{': d += 1
            elif s[j] == '}': d -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

eqs = {b for k in ('problem', 'answer', 'solution') for b in ce_blocks(it[k]) if len(re.findall(r'->|<=>', b)) == 1}
for e in sorted(eqs):
    r = check_equation(e); assert r is None, (e, r); print('式OK:', e)

y = symbols('y', positive=True)
Kp, P = Q(1, 1000), 30
k = sqrt(Q(27, 256) * Kp * P**2); print('k', N(k), 'k^2', N(k**2))
y1 = min(N(r) for r in solve(y - k*(1 - y)**2, y)); print('y1', y1)
y4 = min(N(r) for r in solve(y - k*(Q(9, 10) - y)**2, y)); print('y4', y4)
assert round(float(y1), 3) == 0.198 and round(float(y4), 3) == 0.166
# 問題文の平方根の値で手計算
assert round((1.616 - 1.494)/0.616, 3) == 0.198 and round((1.554 - 1.452)/0.616, 3) == 0.166
assert abs(1.616**2 - 4*0.308**2 - 2.232) < 1e-3 and abs(1.554**2 - 4*0.308*0.249 - 2.109) < 2e-3
# (2) 比 r = xH/xN を変えたときの平衡の y
def y_of_r(r):
    xN = (1 - y)/(1 + r); xH = r*(1 - y)/(1 + r)
    return float(nsolve(y**2 - Kp*P**2*xN*xH**3, y, 0.2))
rs = [1.5, 2, 2.5, 2.9, 3, 3.1, 3.5, 4, 6]
ys = [y_of_r(r) for r in rs]
for r, v in zip(rs, ys): print('r', r, 'y', round(v, 5))
assert max(zip(ys, rs))[1] == 3
x = symbols('x', positive=True); s = symbols('s', positive=True)
sol = solve(diff(x*(s - x)**3, x), x); print(sol); assert s/4 in sol   # x=s は積 0 の端点（最小），最大は x_N2 = s/4, x_H2 = 3s/4
# (5) パージ
w = 0.10/(1 - float(y4)); G = 0.0020/w; loss = (1 - w)*G/0.998
print('w', w, 'G/F', G, 'loss %', loss*100)
assert round(w, 3) == 0.120 and round(loss*100, 1) == 1.5
print('すべて OK')
