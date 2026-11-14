# THE-EQL-116 検証：水相だけで進む平衡 A + B ⇄ P + H2O と分配による生成量の増加
import sympy as sp
x, V, DA, DB, DP, K = sp.symbols('x V D_A D_B D_P K', positive=True)
n0 = sp.Rational(1, 10); Vw = 1
def solve_x(Kv, fA, fB, fP):
    sols = sp.solve(sp.Eq(x*fA*fB, Kv*fP*(n0-x)**2), x)
    sols = [s for s in sols if 0 < s < n0]
    assert len(sols) == 1
    return sols[0]
# (1) 有機溶媒なし
x1 = solve_x(2, 1, 1, 1)
print('(1) x =', sp.N(x1, 5), 'conv =', sp.N(x1/n0*100, 4))
assert round(float(x1/n0*100)) == 15
# 直接の平衡検算
assert abs(float((x1/1)/(((n0-x1)/1)**2)) - 2) < 1e-9
# (2) V = 1.00 L，D_A = D_B = 1.0，D_P = 200
f = lambda D, Vo: Vw + D*Vo
x2 = solve_x(2, f(1,1), f(1,1), f(200,1))
print('(2) x =', sp.N(x2, 5), 'conv =', sp.N(x2/n0*100, 4))
assert round(float(x2/n0*100)) == 73
# 各相の濃度で水相の平衡と分配が成り立つことを確認
Aw = (n0-x2)/f(1,1); Pw = x2/f(200,1)
assert abs(float(Pw/(Aw*Aw)) - 2) < 1e-9
# (3) P の有機相の割合
frac = 200*1/f(200,1)
print('(3) frac =', sp.N(frac*100, 4)); assert round(float(frac*100), 1) == 99.5
# (4) K_app の最大
Kapp = K*(1+DP*V)/((1+DA*V)*(1+DB*V))
dK = sp.diff(sp.log(Kapp.subs({DA:1, DB:1, DP:200, K:2})), V)
Vopt = [s for s in sp.solve(dK, V) if s > 0]
print('(4) Vopt =', Vopt, sp.N(Vopt[0]))
assert len(Vopt) == 1 and abs(float(Vopt[0]) - 0.99) < 1e-12
x4 = solve_x(2, f(1,10), f(1,10), f(200,10))
print('    V=10 L conv =', sp.N(x4/n0*100, 4)); assert round(float(x4/n0*100)) == 58
x4b = solve_x(2, f(1,sp.Rational(99,100)), f(1,sp.Rational(99,100)), f(200,sp.Rational(99,100)))
print('    V=0.99 L conv =', sp.N(x4b/n0*100, 4))
# (5) 少量加えたとき増える条件：d ln Kapp/dV at V=0 = DP - DA - DB
d0 = sp.simplify(sp.diff(sp.log(Kapp), V).subs(V, 0))
print('(5) slope at 0 =', d0); assert sp.simplify(d0 - (DP - DA - DB)) == 0
# 一般の最適値（DA=DB=D）: V = 1/D - 2/DP
D = sp.symbols('D', positive=True)
sol = sp.solve(sp.diff(sp.log(K*(1+DP*V)/(1+D*V)**2), V), V)
print('general Vopt =', sol); assert sp.simplify(sol[0] - (1/D - 2/DP)) == 0
print('OK')
