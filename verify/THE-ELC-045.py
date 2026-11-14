import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R
eqs = ['2H2O -> O2 + 4H+ + 4e-', 'Cu^2+ + 2e- -> Cu', '2H2O + 2e- -> H2 + 2OH-', '4OH- -> O2 + 2H2O + 4e-',
       '2I- -> I2 + 2e-', 'Ag+ + e- -> Ag', 'Cu -> Cu^2+ + 2e-', '2H2 + O2 -> 2H2O']
for e in eqs:
    assert check_equation(e) is None, e
e = R(200, 100)*4825/96500
assert e == R(1, 10)
# 各電極の気体: ('O2', 1/4) は電子 1 mol あたり O2 1/4 mol
cells = {
 'CuSO4': (('O2', R(1,4)), None), 'Na2SO4': (('O2', R(1,4)), ('H2', R(1,2))), 'NaOH': (('O2', R(1,4)), ('H2', R(1,2))),
 'H2SO4': (('O2', R(1,4)), ('H2', R(1,2))), 'KI': (None, ('H2', R(1,2))), 'AgNO3': (('O2', R(1,4)), None),
 'KNO3': (('O2', R(1,4)), ('H2', R(1,2))), 'H2SO4_Cuanode': (None, ('H2', R(1,2))),
}
def residual(a, b):
    g = {'H2': R(0), 'O2': R(0)}
    for c in (a, b):
        for x in cells[c]:
            if x: g[x[0]] += x[1]*e
    h, o = g['H2'], g['O2']
    if h >= 2*o: return ('H2', (h - 2*o)*R(224, 10))
    return ('O2', (o - h/2)*R(224, 10))
res = {k: residual(*v) for k, v in {'ア': ('CuSO4','Na2SO4'), 'イ': ('Na2SO4','NaOH'), 'ウ': ('H2SO4','KI'), 'エ': ('AgNO3','CuSO4'),
                                    'a': ('CuSO4','KNO3'), 'b': ('CuSO4','H2SO4_Cuanode'), 'c': ('CuSO4','CuSO4')}.items()}
for k, v in res.items(): print(k, v[0], v[1], type(v[1]))
for k, want in [('ア', ('O2', R(56,100))), ('イ', ('H2', 0)), ('ウ', ('H2', R(112,100))), ('エ', ('O2', R(112,100)))]:
    print(k, res[k], want); assert res[k][0] == want[0] and res[k][1] == want[1], k
assert res['b'][1] == 0 and res['a'][1] > 0 and res['c'][1] > 0
print('OK')
