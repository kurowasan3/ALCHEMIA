# THE-EQL-050 検証：逐次一次反応 P0->P1->P2->P3（k1, k2=4k1, k3=1.5k1）
import sympy as sp
t = sp.symbols('t', positive=True)
k1 = sp.log(2)/10; k2 = 4*k1; k3 = sp.Rational(3, 2)*k1
x0 = sp.exp(-k1*t)
x1 = k1/(k2-k1)*(sp.exp(-k1*t)-sp.exp(-k2*t))
x2 = k1*k2*(sp.exp(-k1*t)/((k2-k1)*(k3-k1)) + sp.exp(-k2*t)/((k1-k2)*(k3-k2)) + sp.exp(-k3*t)/((k1-k3)*(k2-k3)))
# x2 が微分方程式を満たすか
assert sp.simplify(sp.diff(x2, t) - (k2*x1 - k3*x2)) == 0
assert sp.simplify(sp.diff(x1, t) - (k1*x0 - k2*x1)) == 0
table = {2: (0.002, 0.129, 0.031), 5: (0.024, 0.293, 0.141), 10: (0.119, 0.500, 0.354), 20: (0.399, 0.750, 0.668), 40: (0.800, 0.938, 0.917)}
for T, (val, ala, gly) in table.items():
    a = float((1-x0).subs(t, T)); g = float((1-x0-x1).subs(t, T)); v = float((1-x0-x1-x2).subs(t, T))
    assert (round(v, 3), round(a, 3), round(g, 3)) == (val, ala, gly), (T, v, a, g)
    assert a >= g >= v
# (2)
K1 = 0.693/10; print('k1 =', K1)
# (3) 極大時刻と値
ts = sp.solve(sp.diff(x1, t), t)[0]
assert abs(float(ts) - 20/3) < 1e-9
X0 = float(x0.subs(t, ts)); X1 = float(x1.subs(t, ts))
print('x0(t*)', X0, 'x1max', X1)
assert round(X1, 3) == 0.157 and round(2**(-2/3), 3) == 0.630
K2 = K1*0.630/0.157; print('k2 =', K2); assert f'{K2:.2g}' == '0.28'
print('Gly at t*', 1-0.630-0.157); assert f'{1-0.630-0.157:.2g}' == '0.21'
# (4) k2=100k1
k2b = 100*k1
x1b = k1/(k2b-k1)*(sp.exp(-k1*t)-sp.exp(-k2b*t))
tsb = sp.log(100)/(99*k1)
print('k2=100k1: t*', float(tsb), 'x0', float(x0.subs(t, tsb)), 'x1max', float(x1b.subs(t, tsb)), '1/k2', float(1/k2b))
assert abs(float(x1b.subs(t, tsb)) - 0.0095) < 0.001
print('OK')
