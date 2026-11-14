# THE-EQL-119 検証：A <=> 2B の容器連結と不活性気体
from sympy import Rational as Q, symbols, solve
K = Q(2*10**4); pA0, pB0 = Q(8*10**4), Q(4*10**4)
assert pB0**2/pA0 == K
n0 = pA0 + pB0/2; print('alpha0', float((pB0/2)/n0)); assert (pB0/2)/n0 == Q(1, 5)
y = symbols('y', positive=True)
yb = solve(3*y**2/K + Q(3, 2)*y - n0, y)[0]; pa = yb**2/K
alpha = (yb/2)/(pa + yb/2)
print(float(pa), float(yb), float(pa + yb), float(alpha), float(3*pa/pA0))
assert round(float(pa), -3) == 23000 and round(float(yb), -3) == 21000 and round(float(pa + yb), -3) == 44000 and round(float(alpha), 2) == 0.32
tot3 = pa + yb + Q(10**5)*Q(2, 3); print(float(tot3)); assert round(float(tot3), -4) == 110000
z = solve(y**2/K + y - 6*10**4, y)[0]; pa4 = z**2/K; a4 = (z/2)/(pa4 + z/2)
print(float(pa4), float(z), float(a4))
assert round(float(pa4), -3) == 34000 and round(float(z), -3) == 26000 and round(float(a4), 2) == 0.28
print('OK')
