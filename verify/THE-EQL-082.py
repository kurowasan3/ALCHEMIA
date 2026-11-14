# THE-EQL-082 検証：ビクトル・マイヤー法で PCl5 の解離度と圧平衡定数を求める
from sympy import symbols, solve, Eq, sqrt, N
R, T2 = 8.3e3, 300
Mp = 31.0 + 5*35.5
m = 0.417; n0 = m/Mp
print('n0 =', n0)
assert abs(n0 - 2.000e-3) < 1e-6
def alpha(v, p=1.00e5):     # v [L]
    n = p*v/(R*T2)          # 押し出された空気 = 加熱部で増えた気体粒子
    return n/n0 - 1
def Kp(a, p=1.00e5): return a**2/(1-a**2)*p
a1 = alpha(0.0896); a2 = alpha(0.0707)
K1, K2 = Kp(a1), Kp(a2)
print('(2) 523 K: α =', a1, 'Kp =', K1, ' 見かけの分子量 =', Mp/(1+a1))
print('    473 K: α =', a2, 'Kp =', K2)
assert round(a1,2) == 0.80 and round(a2,2) == 0.42
assert round(K1/1e5,1) == 1.8 and round(K2/1e4,1) == 2.1
assert round(Mp/(1+a1)) == 116
# (3) 高温ほど Kp 大 → 正反応は吸熱
assert K1 > K2
# (4) 試料 2 倍：加熱部の気体は圧力 p のまま → α は同じ，v は 2 倍
a = symbols('a', positive=True)
sol = solve(Eq(a**2/(1-a**2)*1.00e5, K1), a)[0]
print('(4) α（量によらない）=', sol, ' v =', 2*n0*(1+sol)*R*T2/1.00e5*1e3, 'mL')
assert abs(sol - a1) < 1e-9
# 体積一定の容器なら：Kp = (α^2/(1-α)) n0 R T1 / V → n0 が 2 倍で α は小さくなる
V, T1 = 0.100, 523
def a_const_V(n):
    x = symbols('x', positive=True)
    return [s for s in solve(Eq(x**2/(1-x)*n*R*T1/V, K1), x) if 0 < s < 1][0]
print('    定積容器（0.100 L）: n0 で', a_const_V(n0), ' 2n0 で', a_const_V(2*n0))
assert a_const_V(2*n0) < a_const_V(n0)
# (5) p = 0.500e5 Pa
p5 = 0.500e5
a5 = solve(Eq(a**2/(1-a**2)*p5, K1), a)[0]
v5 = n0*(1+a5)*R*T2/p5
print('(5) α =', a5, ' v =', v5*1e3, 'mL')
assert round(float(a5),2) == 0.88 and round(float(v5)*1e3, -1) == 190
print('OK')
