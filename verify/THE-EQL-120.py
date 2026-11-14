# THE-EQL-120 検証：Pd への水素の吸蔵（H2 ⇄ 2H(Pd) の平衡，x = s√p）
import sympy as sp
R, T, V = 8.31e3, 300, 8.31          # 容器1・容器2 とも 8.31 L
c = V/(R*T)                          # 1 Pa あたりの物質量（1 容器）
nPd = 53.0/106
print('Pd', nPd, 'mol,  1容器 1 Pa あたり', c)
# 実験1：容器1 3.50e4 Pa → バルブA を開けて 1.00e4 Pa
n_in = 3.50e4*c; n_gas = 1.00e4*2*c
nH2_abs = n_in - n_gas; x1 = 2*nH2_abs/nPd
print('実験1 吸蔵 H2', nH2_abs, ' x1 =', x1); assert abs(x1-0.200) < 1e-3
# 実験2：A を閉じ（容器2 は 1.00e4 Pa のまま），容器1 を真空にして 8.50e4 Pa → 開けて 4.00e4 Pa
n_in2 = 8.50e4*c + 1.00e4*c; n_gas2 = 4.00e4*2*c
nH2_abs2 = nH2_abs + (n_in2 - n_gas2); x2 = 2*nH2_abs2/nPd
print('実験2 x2 =', x2); assert abs(x2-0.400) < 1e-3
print('x2/x1 =', x2/x1, ' p比 =', 4.00e4/1.00e4, ' sqrt =', (4.00e4/1.00e4)**0.5)
s = x1/(1.00e4)**0.5; K = x1**2/1.00e4
print('s =', s, ' K = x^2/p =', K, K, x2**2/4.00e4)
assert abs(x2**2/4.00e4 - K) < 1e-9
# 実験3：新しい Pd 0.500 mol，両容器真空，容器1 に p0 = 1.25e4 Pa → 開ける
u = sp.symbols('u', positive=True)
sol = sp.solve(sp.Eq(2*c*u**2 + nPd*s*u/2, 1.25e4*c), u)
p3 = float(sol[0]**2); x3 = s*float(sol[0])
abs3 = nPd*x3/2; frac3 = abs3/(1.25e4*c); frac1 = nH2_abs/n_in
print('実験3 p =', p3, ' x =', x3, ' 吸蔵割合', frac3, ' 実験1の割合', frac1)
assert abs(p3-2.50e3) < 1e-6 and abs(frac3-0.60) < 1e-9 and abs(frac1-3/7) < 1e-9
# 比較：もし x = k p（分子のまま溶ける）なら，吸蔵割合は p0 によらない
k = sp.symbols('k', positive=True); P, p0 = sp.symbols('P p0', positive=True)
Ps = sp.solve(sp.Eq(2*c*P + nPd*k*P/2, c*p0), P)[0]
fr = sp.simplify((nPd*k*Ps/2)/(c*p0)); print('x∝p の場合の吸蔵割合', fr); assert p0 not in fr.free_symbols
print('OK')
