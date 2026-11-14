# THE-EQL-122 物質量の変化しない気体平衡に生成物を加えたときの新しい平衡
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, Rational as R, nsimplify
assert check_equation('H2 + I2 <=> 2HI') is None
assert check_equation('2SO2 + O2 <=> 2SO3') is None
# 状態 I: H2, I2 1.00 mol ずつ → HI 1.60 mol
K = R(160, 100)**2 / (R(20, 100) * R(20, 100))
print('K =', K); assert K == 64
# HI 1.00 mol を追加（T, V 一定）。HI が y mol 分解するとして
y = symbols('y')
sols = [s for s in solve((R(26, 10) - 2*y)**2 - K*(R(2, 10) + y)**2, y) if 0 < s < R(13, 10)]
print('y =', sols); y0 = sols[0]
HI = R(26, 10) - 2*y0
print('HI =', HI, 'H2 =', R(2, 10) + y0)
assert HI == R(24, 10) and y0 == R(1, 10)
# 等価な初期量 1.50 mol ずつで比例
assert HI == R(16, 10) * R(15, 10) / 1
xI = R(16, 10) / 2; xII = HI / 3
print('mole fraction HI', xI, xII); assert xI == xII
# 一般式 A(2n+B)/(2n) の確認
A, B, n, z = symbols('A B n z', positive=True)
Kg = A**2 / (n - A/2)**2
cand = A*(2*n + B)/(2*n)
H2 = n - A/2 + (A + B - cand)/2
assert nsimplify(cand**2 / H2**2 - Kg).simplify() == 0
print('all OK')
# (4) SO3 を加えた場合のモル分率の比較（数値例 Kc=10, V=1, 初め SO2 2, O2 1, SO3 を 1 加える）
from sympy import nsolve
x = symbols('x')
def frac(nS, nO, nT, Kc=10):
    # SO3 が 2x 生成（負なら分解）
    sol = [s for s in solve((nT + 2*x)**2 - Kc*(nS - 2*x)**2*(nO - x), x) if s.is_real and -nT/2 < s < nO and 2*s < nS]
    s = sol[0]; tot = nS - 2*s + nO - s + nT + 2*s
    return (nT + 2*s)/tot, s
f1, s1 = frac(2, 1, 0)
nS, nO, nT = 2 - 2*s1, 1 - s1, 2*s1 + 1
f2, _ = frac(nS, nO, nT)
print('SO3 fraction before/after', float(f1), float(f2)); assert f2 > f1
print('SO3 part OK')
