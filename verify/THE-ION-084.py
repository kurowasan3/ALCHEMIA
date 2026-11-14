# THE-ION-084 検証：MX2 型の塩の飽和溶液どうしの混合と析出
import random
from sympy import symbols, solve, Rational as Q, N, nsolve
K = Q(4, 1000)                         # (mol/L)^3
s = symbols('s', positive=True)
s0 = solve(4*s**3 - K, s)[0]; print('pure solubility =', s0); assert s0 == Q(1, 10)
alpha = (Q(1, 10), Q(2, 10)); beta = (Q(4, 10), None); gamma = (None, Q(4, 10))
# beta：[M]=0.40 の飽和溶液 → [X]，gamma：[X]=0.40 の飽和溶液 → [M]
xb = solve(Q(4, 10)*s**2 - K, s)[0]; mg = K/Q(4, 10)**2
beta = (Q(4, 10), xb); gamma = (mg, Q(4, 10)); print('beta', beta, 'gamma', gamma)
assert xb == Q(1, 10) and mg == Q(1, 40)
def mix(a, b, wa=Q(1, 2)): return (wa*a[0] + (1-wa)*b[0], wa*a[1] + (1-wa)*b[1])
for name, (a, b) in {'alpha+beta': (alpha, beta), 'beta+gamma': (beta, gamma), 'alpha+gamma': (alpha, gamma)}.items():
    m, x = mix(a, b); Qv = m*x**2; print(name, 'M', m, 'X', x, 'Q =', float(Qv), '>K' if Qv > K else '<=K')
    assert Qv > K
# (2) 析出量：alpha+beta（等体積）混合液 1 L あたり
m, x = mix(alpha, beta)
t = symbols('t')
ts = nsolve((m - t)*(x - 2*t)**2 - K, t, 0.01); print('(2) precipitated per L =', ts)
assert round(float(ts), 3) == 0.010
# (3) 任意の比で混ぜても Q >= K（数値でランダム検証）
random.seed(1)
for _ in range(20000):
    m1 = random.uniform(0.001, 2); x1 = float((K/m1))**0.5
    m2 = random.uniform(0.001, 2); x2 = float((K/m2))**0.5
    w = random.random()
    mm, xx = w*m1 + (1-w)*m2, w*x1 + (1-w)*x2
    assert mm*xx**2 >= float(K)*(1 - 1e-12)
print('(3) random mixing: always Q >= K')
# (4) alpha を 2 倍に希釈してから beta と等体積混合
a2 = (alpha[0]/2, alpha[1]/2); m, x = mix(a2, beta); Qv = m*x**2
print('(4) diluted alpha + beta: Q =', float(Qv)); assert Qv < K
print('OK')
