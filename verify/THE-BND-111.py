# THE-BND-111 検証：水素化ホウ素の styx（リプスコム）の式の導出と解の列挙
from sympy import symbols, solve, Eq
s, t, y, x, p, q = symbols('s t y x p q')
# 軌道の収支：B は 4 軌道をすべて使う。各 B の基本の B-H で p，追加の末端 H で x，B-H-B は B の軌道 2，BBB は 3，B-B は 2
orb = Eq(p + x + 2*s + 3*t + 2*y, 4*p)
# 電子の収支：価電子 3p + (p+q) = 2 * 結合数（p + x + s + t + y）
ele = Eq(3*p + p + q, 2*(p + x + s + t + y))
hyd = Eq(s + x, q)
sol = solve([orb, ele, hyd], [s, x, y], dict=True)[0]
print('general:', sol)
# Lipscomb の形 s+t=p, t+y=p-q/2, s+x=q に一致するか
assert (sol[s] + t - p).simplify() == 0
assert (sol[y] + t - (p - q/2)).simplify() == 0
def enum(P, Q):
    out = []
    for S in range(0, Q+1):
        for T in range(0, P+1):
            X = Q - S; Yv = P - Q//2 - T
            if Q % 2: continue
            if X < 0 or Yv < 0: continue
            if S + T != P: continue
            if X > P: continue
            out.append((S, T, Yv, X))
    return out
for name, P, Q in [('B2H6', 2, 4), ('B4H10', 4, 6), ('B5H9', 5, 4), ('B6H10', 6, 4), ('B5H11', 5, 6)]:
    print(name, enum(P, Q))
assert enum(2, 4) == [(2, 0, 0, 2)]
assert enum(6, 4) == [(2, 4, 0, 2), (3, 3, 1, 1), (4, 2, 2, 0)]
assert [e for e in enum(6, 4) if e[3] == 0] == [(4, 2, 2, 0)]
assert enum(5, 6) == [(3, 2, 0, 3), (4, 1, 1, 2), (5, 0, 2, 1)]
assert [e for e in enum(5, 6) if e[2] == 0] == [(3, 2, 0, 3)]
# 電子の数の確認：B6H10 (4220)
P, Q, S, T, Yv, X = 6, 4, 4, 2, 2, 0
assert 2*(P + X + S + T + Yv) == 4*P + Q
# かご形 BnHn^2-：骨格電子 = 3n + n + 2 - 2n = 2n + 2
for n, edges in [(6, 12), (12, 30)]:
    print('B%dH%d^2-' % (n, n), 'skeletal electrons', 2*n + 2, 'pairs', n + 1, 'edges', edges)
# 二重項の数：BnHn^2- の総価電子
