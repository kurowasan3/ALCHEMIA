# THE-GAS-074 verification: protein molar mass from 2D gaseous film
from sympy import symbols, solve, Rational as R, N
I, B = symbols('I B')
S1, p1, S2, p2 = R(200, 1000), R(100, 1000) * R(1, 1000), R(100, 1000), R(300, 1000) * R(1, 1000)
sol = solve([I + B * p1 - p1 * S1, I + B * p2 - p2 * S2], [I, B])
print(sol); assert sol[I] == R(3, 200000) and sol[B] == R(1, 20)
RT = 8.31 * 298; m = 0.110e-3; NA = 6.0e23
n = float(sol[I]) / RT; M = m / n; a0 = float(sol[B]) / (n * NA)
print('RT', RT, 'n', n, 'M', M, 'a0 nm2', a0 * 1e18)
assert round(M, -3) == 18000 and round(a0 * 1e18) == 14
for S, p in ((0.200, 1.00e-4), (0.100, 3.00e-4)):
    nn = p * S / RT; print('naive M', m / nn)
    assert abs(p - float(sol[I]) / (S - float(sol[B]))) < 1e-12
assert round(m / (2.00e-5 / RT), -3) == 14000 and round(m / (3.00e-5 / RT), -2) == 9100
print('pi at 0.060', float(sol[I]) / (0.060 - float(sol[B])))
assert abs(float(sol[I]) / (0.060 - 0.050) - 1.5e-3) < 1e-9
print('OK')
