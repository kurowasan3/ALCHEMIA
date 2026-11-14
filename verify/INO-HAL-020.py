# 検証スクリプト（グループ g11）  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/INO-HAL-020/item.json'
import sys, json, re, math
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
assert check_equation('2XeF2 + 2H2O -> 2Xe + 3HF + O2') is not None  # 検査器が誤りを検出できること
def sig(x, n):
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
import sympy as sp
# XeF4 の加水分解：Xe の行き先が未定なら 1 パラメータの不定性がある
x, y, z, w, u = sp.symbols('x y z w u')  # 6XeF4 + x H2O -> y Xe + z XeO3 + w HF + u O2
# 電子：Xe(+4)->0 で 4y 受け取り，+4->+6 で 2z 出し，O で 4u 出す
sol = sp.solve([sp.Eq(6, y+z), sp.Eq(w, 24), sp.Eq(2*x, w), sp.Eq(x, 3*z+2*u), sp.Eq(4*y, 2*z+4*u)], [x, y, w, u], dict=True)[0]
print('一般解', sol)  # z が自由
assert sol[u].free_symbols == {z}  # O2 の係数が XeO3 の係数 z に依存して決まらない
s4 = {k: v.subs(z, 2) for k, v in sol.items()}
print('z=2', s4); assert s4 == {x: 12, y: 4, w: 24, u: 3}
s_other = {k: v.subs(z, 1) for k, v in sol.items()}; print('z=1 でもつり合う', s_other)
assert check_equation('6XeF4 + 12H2O -> 4Xe + 2XeO3 + 24HF + 3O2') is None
assert check_equation('12XeF4 + 24H2O -> 10Xe + 2XeO3 + 48HF + 9O2') is None  # 別の割合でもつり合う（不定性の確認）
# 酸化数
for f, n in (('XeF2', 2), ('XeF4', 4), ('XeF6', 6)): assert n == {'XeF2': 2, 'XeF4': 4, 'XeF6': 6}[f]
# (5)
M = 131 + 2*19.0; n = 3.38/M; print(M, n); assert abs(n - 0.0200) < 1e-4
gas = n*1.5*22.4*1000; print('気体 mL', gas); assert sig(gas, 3) == 672
naoh = 2*n/0.500*1000; print('NaOH mL', naoh); assert sig(naoh, 3) == 80.0
# (6)
g2 = sp.Rational(3, 2); g4 = sp.Rational(4, 6) + sp.Rational(3, 6); g6 = 0
r = [g2*6, g4*6, g6]; print('比', r); assert r == [9, 7, 0]
assert sp.Rational(1, 2)/g2 == sp.Rational(1, 3) and sp.Rational(1, 2)/g4 == sp.Rational(3, 7)
print('OK')
