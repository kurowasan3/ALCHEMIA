# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem2/outh/items/THE-EQL-055/item.json'
it = json.load(open(ITEM, encoding='utf-8'))


def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0:
            return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i + 4:j - 1]); i = j


eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)
# 小問数と答えの番号の一致
nq = it['problem'].count('\\item')
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na == 0), (nq, na)
print('小問', nq, '答え', na)


def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
import sympy as sp
R = sp.Rational
x, y = sp.symbols('x y')


def ester(a, b, w):
    s = [t for t in sp.solve(sp.Eq(x * (x + w), 4 * (a - x) * (b - x)), x) if 0 < t < a]
    assert len(s) == 1
    return float(s[0])


x0 = ester(1, 2, 0)
print('副反応なし', x0); assert sig(x0, 3) == 0.845
ys = sp.solve(sp.Eq(R(3, 4) * (R(3, 4) + y), 4 * R(1, 4) * (2 - R(3, 4) - 2 * y)), y)
print('y =', ys); assert ys == [R(1, 4)]
yv = ys[0]
K = R(3, 4) * (R(3, 4) + yv) / (R(1, 4) * (2 - R(3, 4) - 2 * yv))
assert K == 4
M = 4 * 12 + 8 * 1.0 + 2 * 16
assert M == 88 and sig(0.750 * M, 3) == 66.0
xi = ester(1, R(3, 2), 0); xii = ester(1, 2, R(1, 4))
print('(i)', xi, '(ii)', xii)
assert sig(xi, 3) == 0.785 and sig(xii, 3) == 0.816
assert abs(ester(1, R(3, 2), R(1, 4)) - 0.75) < 1e-12
assert (x0 - xi) > (x0 - xii)
print('減少', x0 - xi, x0 - xii, '和', (x0 - xi) + (x0 - xii), '実際', x0 - 0.75)
assert check_equation('CH3COOH + C2H5OH <=> CH3COOC2H5 + H2O') is None
print('THE-EQL-055 検証OK')
