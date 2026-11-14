# THE-EQL-099 検証：生成物が中間体を2倍にする表面反応（連鎖分岐型の自己触媒）
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
k3, X = 5.0e6, 5.0e-10
kp = k3*X
print("k' =", kp); assert abs(kp - 2.5e-3) < 1e-15
t1000 = 3*2.30/kp; t2 = 0.69/kp
print('t(1000x)=', t1000, ' t(2x)=', t2); assert f"{t1000:.1e}" == '2.8e+03' and f"{t2:.1e}" == '2.8e+02'
# 完全な連立微分方程式で確認（(2) が速い：k2[Y] = 2 /s ≫ k3[X]）
k2, Y = 1.0e9, 2.0e-9
def f(t, z):
    XY, P = z
    v2 = k2*XY*Y; v3 = k3*X*P
    return [2*v3 - v2, v2 - v3]
P0 = 1e-15
sol = solve_ivp(f, [0, 4000], [0, P0], rtol=1e-10, atol=1e-30, dense_output=True)
tt = np.linspace(100, 4000, 50)
ratio = sol.sol(tt)[1]/P0
growth = np.log(ratio)/tt
print('effective rate (late)', growth[-1]); assert abs(growth[-1] - kp)/kp < 0.01
# 1000 倍になる時間
from scipy.optimize import brentq
tx = brentq(lambda t: sol.sol(t)[1]/P0 - 1000, 100, 4000); print('numerical t1000', tx)
assert abs(tx - 2763)/2763 < 0.01
# 分岐しない場合 (X + XY2 -> XY + Q)：XY2 の量は一定
def g(t, z):
    XY, P = z
    v2 = k2*XY*Y; v3 = k3*X*P
    return [v3 - v2, v2 - v3]
s2 = solve_ivp(g, [0, 4000], [0, P0], rtol=1e-10, atol=1e-30)
tot = s2.y[0][-1] + s2.y[1][-1]
print('non-branching total', tot/P0); assert abs(tot/P0 - 1) < 1e-6
# Q の生成速度 = k3[X][XY2]_0
print('Q rate', kp*P0)
