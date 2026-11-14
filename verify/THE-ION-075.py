from sympy import Rational as R, log, N, solve, symbols
KIn = R(4, 10)*R(1, 10**9)
cv = R(2, 10**6)
def H_at(C): return KIn*(C - cv)/cv
C1, C2 = R(2, 10**5), R(2, 10**4)
pH1 = -log(H_at(C1), 10); pH2 = -log(H_at(C2), 10)
print('pH1', N(pH1, 5), 'pH2', N(pH2, 5))
assert round(float(pH1), 1) == 8.4 and round(float(pH2), 1) == 7.4
Ka = R(1, 10**6)
f = lambda H: H/(H + Ka)
f1, f2 = f(H_at(C1)), f(H_at(C2))
print('残存率', N(f1*100, 4), N(f2*100, 4))
assert round(float(f1*100), 2) == 0.36 and round(float(f2*100), 1) == 3.8
# 近似の確認
nHB = R(1, 10)*R(20, 1000)
V = R(40, 1000)
In_consumed = cv*V
OH1 = R(1, 10**14)/H_at(C1)*V
print('指示薬の消費/HB', N(In_consumed/nHB, 3), 'OH-/HB', N(OH1/nHB, 3))
assert In_consumed/nHB < R(1, 10**4) and OH1/nHB < R(1, 10**4)
C = symbols('C', positive=True)
Hmax = R(2, 98)*Ka
Cmax = solve(KIn*(C - cv)/cv - Hmax, C)[0]
print('Cmax', N(Cmax, 5))
assert abs(float(Cmax) - 1.04e-4) < 1e-6
print('OK')
