"""INO-NIT-016 検証：金属と硝酸の反応で酸化剤としてはたらく硝酸の割合と合金の組成（SymPy／chemlib）"""
import sys
from sympy import symbols, Rational as R, solve, linsolve
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

eqs = {'3Cu + 8HNO3 -> 3Cu(NO3)2 + 2NO + 4H2O': (8, 2), 'Cu + 4HNO3 -> Cu(NO3)2 + 2NO2 + 2H2O': (4, 2),
       '3Ag + 4HNO3 -> 3AgNO3 + NO + 2H2O': (4, 1), 'Ag + 2HNO3 -> AgNO3 + NO2 + H2O': (2, 1)}
for e, (tot, ox) in eqs.items():
    r = check_equation(e); print(e, r or 'OK', '酸化剤の割合', R(ox, tot)); assert r is None
    assert R(ox, tot) in (R(1, 4), R(1, 2))
# 一般の n 価の金属 M：3M + 4n HNO3 -> 3M(NO3)n + n NO + 2n H2O
n = symbols('n', positive=True)
print('希硝酸 一般', (n) / (4 * n), '濃硝酸 一般', n / (2 * n))
assert (n) / (4 * n) == R(1, 4) and n / (2 * n) == R(1, 2)
# 一般式 3M + 4n HNO3 -> 3M(NO3)n + n NO + 2n H2O を具体の価数で確認（n=1 Ag, 2 Cu, 3 Fe）
for M, k in (('Ag', 1), ('Cu', 2), ('Fe', 3)):
    salt = f'{M}NO3' if k == 1 else f'{M}(NO3){k}'
    for e in (f'3{M} + {4*k}HNO3 -> 3{salt} + {k}NO + {2*k}H2O', f'{M} + {2*k}HNO3 -> {salt} + {k}NO2 + {k}H2O'):
        r = check_equation(e.replace(' 1NO', ' NO').replace(' 1NO2', ' NO2').replace(' 1H2O', ' H2O')); print(e, r or 'OK'); assert r is None
# 合金：x mol Cu, y mol Ag
x, y = symbols('x y')
sol = solve([R(635, 10) * x + 108 * y - R(686, 100), 2 * x + y - 3 * R(896, 1000) / R(224, 10)], [x, y])
print(sol); assert sol[x] == R(4, 100) and sol[y] == R(4, 100)
mCu, mAg = sol[x] * R(635, 10), sol[y] * 108; print('Cu g', float(mCu), 'Ag g', float(mAg))
nNO = R(896, 1000) / R(224, 10)
HNO3_dil = 4 * nNO; print('希硝酸の消費', HNO3_dil); assert HNO3_dil == R(16, 100)
# 組成から直接数えても一致
assert sol[x] * R(8, 3) + sol[y] * R(4, 3) == HNO3_dil
ne = 2 * sol[x] + sol[y]
NO2 = ne; HNO3_conc = 2 * ne
print('濃硝酸：NO2 mol', NO2, 'L', float(NO2 * R(224, 10)), 'HNO3', HNO3_conc)
assert round(float(NO2 * R(224, 10)), 2) == 2.69 and HNO3_conc == R(24, 100)
assert sol[x] * 4 + sol[y] * 2 == HNO3_conc
print('INO-NIT-016 all OK')
