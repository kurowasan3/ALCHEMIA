# 検証スクリプト THE-SOL-042  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-SOL-042/item.json'
import sys, json, re
from math import log10, floor
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
def check_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        print('式', e, '->', 'OK' if r is None else r)
        assert r is None, (e, r)
    return sorted(eqs)
def sig(x, n):
    return float(f'{x:.{n}g}')
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
def close(a, b, rel=1e-9):
    assert abs(a - b) <= rel * max(abs(a), abs(b)), (a, b)
check_structure()
from sympy import Rational as Q, symbols, solve, Eq
eqs = check_equations(); assert 'CuSO4 + 5H2O -> CuSO4.5H2O' in eqs
f = Q(160, 250); w = Q(90, 250); S = Q(25, 100)
c0, w0 = Q(15), Q(85)
m, x, y = symbols('m x y')
m1 = solve(Eq(c0 + m, S*w0), m)[0]; y1 = solve(Eq((c0 + f*y)/(w0 + w*y), S), y)[0]
print('(1) 無水物', m1, float(m1), '五水和物', float(y1)); assert m1 == Q(25, 4) and sig(float(y1), 3) == 11.4
xs = solve(Eq((c0 + m - f*x)/(w0 - w*x), S), x)[0]; print('(2) x(m) =', xs)
x20 = xs.subs(m, 20); print('(2) x', x20, float(x20)); assert x20 == 25
cu_cr = f*x20; print('結晶中 CuSO4', cu_cr, '超過分', 20 - m1); assert cu_cr == 16 and cu_cr > 20 - m1
rem_c, rem_w = c0 + 20 - cu_cr, w0 - w*x20; print('残液', rem_c, rem_w); assert rem_c/rem_w == S
xmax = w0/w; mmax = solve(Eq(xs, xmax), m)[0]; print('(3) x_max', float(xmax), 'm', float(mmax))
assert sig(float(mmax), 3) == 136 and sig(float(xmax), 3) == 236 and f*xmax == c0 + mmax
print('(4) 傾き', float(1/(f - w*S)), '五水和物の場合の傾き 1')
assert sig(float(1/(f - w*S)), 3) == 1.82
# 五水和物を飽和をこえて加えた場合：溶液は変化せず溶け残りは y - y1
# 数値で平衡を総当たり確認：無水物 m を加えたとき，溶液の濃度が S を超えない最小の析出量
for mv in (0, 5, 6.25, 10, 20, 100, 136, 140):
    xv = max(0.0, float(xs.subs(m, mv))); xv = min(xv, float(xmax))
    cw = float(w0 - w*xv); cc = float(c0 + mv - f*xv)
    if cw > 1e-9: assert cc/cw <= 0.25 + 1e-9
    print('  m', mv, 'x', round(xv, 2), '残液の水', round(cw, 2))
print('すべて一致')
