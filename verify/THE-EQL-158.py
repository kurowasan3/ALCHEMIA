# THE-EQL-158 解離吸着のラングミュア型モデルとデータからの判定
from sympy import symbols, solve, sqrt, simplify, Eq, Rational as Rr
import numpy as np
th, K, P = symbols('theta K P', positive=True)
# (1) ka P (1-θ)^2 = kd θ^2 → θ
sol = solve(Eq(K*P*(1-th)**2, th**2), th); print('theta solutions:', [simplify(s) for s in sol])
# データ生成（真値 N0 = 2.0e-4 mol サイト，K = 4.0e-4 /Pa）
N0, Kv = 2.0e-4, 4.0e-4
Ps = np.array([25., 100., 400., 1600.])
theta = np.sqrt(Kv*Ps)/(1+np.sqrt(Kv*Ps))
n = N0*theta/2          # 吸着した H2 の物質量
print('n (mol) =', n)
# 分子状モデル：P/n vs P
y1 = Ps/n; print('P/n =', y1, ' slopes', np.diff(y1)/np.diff(Ps))
# 解離モデル：√P/n vs √P
x2 = np.sqrt(Ps); y2 = x2/n; print('sqrtP/n =', y2, ' slopes', np.diff(y2)/np.diff(x2))
slope = np.diff(y2)[0]/np.diff(x2)[0]; inter = y2[0]-slope*x2[0]
print('slope', slope, '-> N0 =', 2/slope, ' intercept', inter, '-> K =', (slope/inter)**2)
# (3) θ=1/2 となる圧力，低圧で P を 1/4 にしたときの吸着量
print('P(θ=0.5) =', 1/Kv)
lowP = 1e-2
r = (np.sqrt(Kv*lowP/4)/(1+np.sqrt(Kv*lowP/4)))/(np.sqrt(Kv*lowP)/(1+np.sqrt(Kv*lowP))); print('low P, P/4 -> ratio', r)
# 表面積：サイト 1 個 6.0e-20 m^2
print('area m^2 =', N0*6.02e23*6.0e-20)
