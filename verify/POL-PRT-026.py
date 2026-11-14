# 検証スクリプト POL-PRT-026  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/POL-PRT-026/item.json'
import json, re
from math import log10
from sympy import symbols, simplify, solve, nsolve
it = json.load(open(ITEM, encoding='utf-8'))
n_items = len(re.findall(r'\\item\b', it['problem'])); n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans == 4, (n_items, n_ans)
eqs = [b for k in ('problem', 'answer', 'solution') for b in re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', it[k])
       if len(re.findall(r'->|<=>|<-', b)) == 1]
assert eqs == []
def sig(x, n):
    from math import floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
# 記号的な関係
C, h, kZ, kN, kZp = symbols('C h kZ kN kZp', positive=True)
Z = kZ*C/h; N = kN*C/h; A = kZp*Z/h; kNp = A*h/N
assert simplify((Z + N)*h/C - (kZ + kN)) == 0
assert simplify(kZ*kZp - kN*kNp) == 0
K2 = A*h/(Z + N); assert simplify(1/K2 - (1/kZp + 1/kNp)) == 0
# 数値
K1, K2v, kNv = 10**-2.34, 10**-9.69, 10**-7.80
kZv = K1 - kNv; r = kZv/kNv; print('Z/N', r); assert sig(r, 2) == 2.9e5
# k'Z, k'N を厳密に：1/K2 = 1/kZp + 1/kNp, kZ kZp = kN kNp → kNp = kZ kZp/kN
kzp = K2v*(1 + kNv/kZv)  # 1/K2 = (1/kzp)(1 + kN/kZ)
knp = kZv*kzp/kNv; print("k'Z", -log10(kzp), "k'N", -log10(knp)); assert round(-log10(knp), 2) == 4.23
# 近似値
assert round(2.34 + 9.69 - 7.80, 2) == 4.23
# 比べる倍率
print('倍率', kZv/knp, kNv/kzp); assert abs(log10(kZv/knp) - 1.89) < 0.01 and abs(log10(kNv/kzp) - 1.89) < 0.01
# 分布を微視的定数から計算して見かけの定数を再現
def dist(hv):
    c = 1.0; z = kZv*c/hv; n = kNv*c/hv; a = kzp*z/hv
    tot = c + z + n + a
    return c/tot, z/tot, n/tot, a/tot
for pH in (2.0, 6.0, 10.0):
    hv = 10**-pH; c_, z_, n_, a_ = dist(hv)
    assert abs((z_ + n_)*hv/c_ - K1)/K1 < 1e-9 and abs(a_*hv/(z_ + n_) - K2v)/K2v < 1e-9
    assert abs(z_/n_ - r)/r < 1e-9
pI = (2.34 + 9.69)/2; print('pI', pI); assert round(pI + 1e-9, 2) == 6.02
c_, z_, n_, a_ = dist(10**-pI); print('C', c_, 'N', n_, 'A', a_)
assert abs(c_ - a_) / c_ < 1e-6 and sig(c_, 1) == 2e-4 and sig(n_, 1) == 3e-6
assert abs(10**0.32 - 2.1) < 0.02 and abs(10**0.46 - 2.9) < 0.02
print('全検証OK:', ITEM)
