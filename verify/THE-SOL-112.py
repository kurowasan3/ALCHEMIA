# THE-SOL-112 検証：NaCl-KCl-水（25℃）の溶解度図による析出量・最小の水・逆問題
from sympy import Rational as R, symbols, solve, nsimplify
# 溶解度図（水 100 g あたり g）：KCl 飽和線 K = 36 - 2N/3 (0<=N<=30)，NaCl 飽和線 N = 36 - 3K/8 (0<=K<=16)，共飽和点 B(30,16)
def Ksat(N): return R(36) - R(2,3)*N
def Nsat(K): return R(36) - R(3,8)*K
def equil(a, b, W):
    """NaCl a g, KCl b g, 水 W g -> (残る NaCl, 残る KCl)"""
    N, K = R(100)*a/W, R(100)*b/W
    # 領域 III
    if N > 30 and K > 16:
        return (a - R(30)*W/100, b - R(16)*W/100)
    # KCl のみ析出（溶液は NaCl をすべて含む）
    if N <= 30 and K > Ksat(N):
        return (0, b - Ksat(N)*W/100)
    if K <= 16 and N > Nsat(K):
        return (a - Nsat(K)*W/100, 0)
    return (0, 0)
# (1)
r = equil(R(20), R(30), R(50)); print('(1)', r); assert r == (5, 22)
# (2) 最小の水
W = symbols('W', positive=True)
sol = solve(R(100)*30/W - Ksat(R(100)*20/W), W); print('(2) W =', sol, [float(s) for s in sol])
Wmin = sol[0]
assert abs(float(Wmin) - 120.4) < 0.1
# ray slope check: K/N = 1.5 > 16/30 so it meets KCl line
print('  separately: KCl', 30/0.36, 'NaCl', 20/0.36)
# 120 g では溶け切り，119 g では析出
print('  W=121 ->', equil(R(20), R(30), R(121)), ' W=119 ->', [float(x) for x in equil(R(20), R(30), R(119))])
# (3) 逆問題：合計 60 g，水 100 g，残り 18.0 g
sols = []
for i in range(0, 60001):
    a = R(i, 1000); b = 60 - a
    rN, rK = equil(a, b, R(100))
    if rN + rK == 18:
        sols.append((a, b, rN, rK))
print('(3) solutions:', [(float(a), float(b), float(x), float(y)) for a, b, x, y in sols])
assert [(a, b) for a, b, _, _ in sols] == [(R(18), R(42)), (R(252, 5), R(48, 5))]
# 残り 14 g のときは NaCl 30〜44 g の全範囲
cnt = [float(a) for a in (R(i, 10) for i in range(0, 601)) if sum(equil(a, 60 - a, R(100))) == 14]
print('  residue 14 g: NaCl from', min(cnt), 'to', max(cnt), 'count', len(cnt))
