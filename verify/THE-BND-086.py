# THE-BND-086 検証：閃亜鉛鉱型の充填率と半径比（陽イオン-陰イオン接触と陰イオンどうしの接触の切りかわり）
import math
from sympy import symbols, solve, sqrt, Rational as Q, nsimplify
def phi(r, pi=math.pi, s2=math.sqrt(2), s3=math.sqrt(3)):
    rc = math.sqrt(6) / 2 - 1
    if r >= rc:
        a = 4 * (1 + r) / s3
    else:
        a = 2 * s2
    return 4 * (4 / 3) * pi * (1 + r**3) / a**3
r = symbols('r', positive=True)
rc = solve(sqrt(2) * 4 * (1 + r) / sqrt(3) / 2 - 2, r)[0]
print('critical', rc, float(rc))
assert abs(float(rc) - 0.2247) < 1e-3
for rr in (0.40, 0.15, float(rc)):
    exact = phi(rr); approx = phi(rr, 3.14, 1.414, 1.732)
    print(rr, exact, approx)
p40, p15, pc = phi(0.40), phi(0.15), phi(float(rc))
assert f'{p40:.2g}' == '0.53' and f'{phi(0.40,3.14,1.414,1.732):.2g}' == '0.53'
assert f'{p15:.2g}' == '0.74' and f'{phi(0.15,3.14,1.414,1.732):.2g}' == '0.74'
assert f'{pc:.2g}' == '0.75' and f'{phi(float(rc),3.14,1.414,1.732):.2g}' == '0.75'
# 単調性：r<rc で増加，r>rc で減少 → rc で最大
import numpy as np
xs = np.linspace(0.05, 0.6, 2000)
ys = [phi(x) for x in xs]
imax = int(np.argmax(ys)); print('max at', xs[imax], ys[imax])
assert abs(xs[imax] - float(rc)) < 1e-3
# 最近接の陰イオンどうしの距離（r=0.40）
a = 4 * 1.40 / math.sqrt(3)
print('X-X distance / r-', a / math.sqrt(2), '> 2')
assert a / math.sqrt(2) > 2
print('OK')
