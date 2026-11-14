# 検証スクリプト（グループ g03）  python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-EQL-151/item.json'
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

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
import sympy as sp
# データの作り方の検証：仕込み 1:2，K=4 の真の平衡から積分値を作る
x = sp.symbols('x', positive=True)
sol = [s for s in sp.solve(sp.Eq(x**2, 4*(1-x)*(2-x)), x) if s < 1][0]
xe = float(sol); e = 2 - xe; c = 1 - xe
scale = 30.0 / (3*(c + xe))
I = dict(o412=2*xe*scale, o369=2*e*scale, ac=3*(c+xe)*scale, me=3*(e+xe)*scale)
print({k: round(v, 3) for k, v in I.items()})
assert round(I['o412'], 1) == 16.9 and round(I['o369'], 1) == 23.1 and round(I['me'], 1) == 60.0
# 生徒の解き方（与えた積分値から）
est = 16.9/2; eth = 23.1/2; acid = 30.0/3 - est
assert abs(est-8.45) < 1e-9 and abs(acid-1.55) < 1e-9 and sig(eth, 3) == 11.6
assert abs(3*(eth+est) - 60.0) < 1e-9          # (2) 確認
assert abs((acid+est) - 10.0) < 1e-9 and abs((eth+est) - 20.0) < 1e-9   # (3) 1:2
K = est*est/(acid*eth)
print('K', K); assert sig(K, 2) == 4.0
# (5)
y = sp.symbols('y')
r = [s for s in sp.solve(sp.Eq(y*(1+y), 4*(1-y)**2), y) if 0 < s < 1]
yv = float(r[0]); print('x', yv, 'ratio', yv/(1-yv))
assert abs(yv - (9-33**0.5)/6) < 1e-12 and round(yv, 3) == 0.543 and sig(yv/(1-yv), 2) == 1.2
# 水なし 1:1 なら x=2/3, 比 2.0
z = [s for s in sp.solve(sp.Eq(y**2, 4*(1-y)**2), y) if 0 < s < 1][0]
assert z == sp.Rational(2, 3)
print('OK')
