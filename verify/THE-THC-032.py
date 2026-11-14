# 検証スクリプト（グループ g03）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, factor, simplify, solve, diff, sqrt, Rational as Q, N, Eq
ITEM = '/home/claude/chem4/outh/items/THE-THC-032/item.json'
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
            if '\\text' in b or '\\mathrm' in b:
                print('一般式（検査対象外）:', b); continue
            r = check_equation(b); assert r is None, (b, r); print('式OK:', b)

def sig(x, n=2):
    from math import log10, floor
    x = float(x); return round(x, -int(floor(log10(abs(x)))) + (n - 1))

r, s, k, x = symbols('r s k x', positive=True)
g4 = k / (r + s) - k / 4 * (1 / r + 1 / s)
assert simplify(g4 + k * (r - s)**2 / (4 * r * s * (r + s))) == 0          # (1)
g5 = k / (r + s) - k / 5 * (1 / r + 1 / s)
assert simplify(g5 + k * (r**2 - 3 * r * s + s**2) / (5 * r * s * (r + s))) == 0
roots = sorted(solve(Eq(x**2 - 3 * x + 1, 0), x), key=lambda v: float(v))
assert [round(float(v), 2) for v in roots] == [0.38, 2.62]                   # (2) 範囲
crit = solve(diff(g5.subs(r, x * s), x), x)
assert len(crit) == 1 and simplify(crit[0] - (1 + sqrt(5)) / 4) == 0 and round(float(crit[0]), 2) == 0.81
d2 = diff(g5.subs(r, x * s), x, 2).subs(x, crit[0]).subs(k, 1).subs(s, 1)
assert float(d2) < 0                                                          # 最大
K = Q(240)
G = lambda a, b: K / (a + b) - K / 5 * (1 / a + 1 / b)
rad = {'Li': Q('0.076'), 'Na': Q('0.102'), 'K': Q('0.138'), 'F': Q('0.133'), 'I': Q('0.220'), 'ClO4': Q('0.240')}
na, kk = G(rad['Na'], rad['ClO4']), G(rad['K'], rad['ClO4'])
assert sig(na) == 31 and sig(kk) == 87 and kk > na                            # (3)
lif, lii = G(rad['Li'], rad['F']), G(rad['Li'], rad['I'])
assert lif > 0 and lii < 0 and sig(lif) == 160 and sig(lii) == -39            # (4)
assert round(float(rad['Li'] / rad['F']), 2) == 0.57 and round(float(rad['Li'] / rad['I']), 2) == 0.35
print('NaClO4', N(na, 4), 'KClO4', N(kk, 4), 'LiF', N(lif, 4), 'LiI', N(lii, 4))
print('THE-THC-032 OK')
