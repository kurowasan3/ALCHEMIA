# THE-ION-104 検証：共通イオン効果で近似が使えない場合（2次方程式）と近似の誤差
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve, sqrt, Rational as R, N
K = R(3, 10)*R(10)**-9      # 3.0e-10
s = symbols('s', positive=True)
s0 = sqrt(K); print('s0 =', N(s0))
assert round(float(s0)/1e-5, 1) == 1.7
c = R(2)*R(10)**-5
x = solve(s*(s + c) - K, s)[0]; print('exact [Ag+] =', N(x))
assert abs(float(x) - 1.0e-5) < 1e-12
approx = K/c; print('approx =', N(approx), 'rel.err =', N((approx - x)/x))
assert abs(float(approx) - 1.5e-5) < 1e-12
assert float(approx) < float(s0)   # 近似値は s0 より小さいが誤差 50%
# 相対誤差 = x/c を一般に確認
cc = symbols('cc', positive=True)
xx = (-cc + sqrt(cc**2 + 4*K))/2
assert abs(float(((K/cc - xx)/xx - xx/cc).subs(cc, R(37, 10**6)))) < 1e-12
# 誤差 10% 以下：x <= 0.1c  <=> K <= 0.11 c^2
cmin = sqrt(K/R(11, 100)); print('c_min =', N(cmin))
assert round(float(cmin)/1e-5, 1) == 5.2
x_at = (-cmin + sqrt(cmin**2 + 4*K))/2
assert abs(float(x_at/cmin) - 0.1) < 1e-9
print('sqrt(27.3)=', N(sqrt(R(273, 10))))
assert check_equation('AgCl <=> Ag+ + Cl-') is None
print('OK')
