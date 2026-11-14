# THE-GAS-027 検証：空気＋揮発性物質の等温圧縮データから飽和蒸気圧と物質量を逆算
from sympy import Rational as Q, symbols, solve, N
R = Q(83,10)*1000; T = 300; RT = R*T                  # 2.49e6 Pa L/mol
na, nx, ps = Q(2,100), Q(5,100), 6*10**4
def P(V):
    if nx*RT/V >= ps:           # 液体が存在
        return na*RT/V + ps, True
    return (na+nx)*RT/V, False
data = {}
for V in [Q(1), Q(3,2), Q(3), Q(5)]:
    p, liq = P(V); data[V] = p; print(V, N(p, 6), liq, 'PV=', N(p*V, 6))
given = {Q(1): Q(1098,1000)*10**5, Q(3,2): Q(9320,10000)*10**5, Q(3): Q(5810,10000)*10**5, Q(5): Q(3486,10000)*10**5}
for V in given: assert abs(given[V] - data[V]) / data[V] < 6e-4, V
# 学生の手順：PV の値から液体の有無を判定
PV = {V: given[V]*V for V in given}; print({k: N(v) for k, v in PV.items()})
# 1.00, 1.50 L：液体あり → P = A/V + ps
A, s = symbols('A s')
sol = solve([A/1 + s - given[Q(1)], A/Q(3,2) + s - given[Q(3,2)]], [A, s]); print(sol)
assert abs(sol[s] - 6.0e4) < 1 and abs(sol[A] - 4.98e4) < 1
n_air = sol[A]/RT; n_tot = (given[Q(3)]*3)/RT; print('n_air', N(n_air), 'n_tot', N(n_tot), N((given[Q(5)]*5)/RT))
assert abs(n_air - 0.020) < 1e-5 and abs(n_tot - 0.070) < 2e-4
Vstar = (n_tot - n_air)*RT/sol[s]; Pstar = n_tot*RT/Vstar
print('V*', N(Vstar), 'P*', N(Pstar)); assert abs(Vstar - 2.08) < 0.01 and abs(Pstar - 8.4e4) < 0.05e4
# 2.00 L での液体
nliq = Q(5,100) - sol[s]*2/RT; print('2.00 L liquid', N(nliq)); assert abs(nliq - 1.8e-3) < 0.05e-3
# 2.00 L での全圧
print('P(2.00)', N(n_air*RT/2 + sol[s]))
print('OK')
