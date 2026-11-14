# THE-BND-102 verification: MO filling, bond orders, unpaired electrons
from fractions import Fraction as Fr
def fill(levels, n):
    """levels: list of (kind, degeneracy) low->high. Returns bonding, antibonding, unpaired."""
    b = a = unp = 0
    for kind, g in levels:
        k = min(n, 2*g); n -= k
        if kind == 'b': b += k
        else: a += k
        unp += k if k <= g else 2*g - k
    assert n == 0
    return b, a, unp
H = [('b',1),('a',1)]
res = {}
for name, ne in [('H2+',1),('H2',2),('H2-',3),('He2+',3),('He2',4)]:
    b,a,u = fill(H, ne); res[name] = Fr(b-a,2); print(name, b, a, res[name])
assert res == {'H2+':Fr(1,2),'H2':1,'H2-':Fr(1,2),'He2+':Fr(1,2),'He2':0}
O = [('b',1),('a',1),('b',1),('b',2),('a',2),('a',1)]
bo = {}; up = {}
for name, ne in [('O2',12),('O2+',11),('O2-',13),('O2^2-',14)]:
    b,a,u = fill(O, ne); bo[name] = Fr(b-a,2); up[name] = u; print(name, b, a, bo[name], 'unpaired', u)
assert bo == {'O2':2,'O2+':Fr(5,2),'O2-':Fr(3,2),'O2^2-':1}
assert up == {'O2':2,'O2+':1,'O2-':1,'O2^2-':0}
order = sorted(bo, key=lambda k: -bo[k]); print('shortest first', order)
assert order == ['O2+','O2','O2-','O2^2-']
# literature: H2 436 kJ/mol 74 pm; H2+ ~2.65 eV = 256 kJ/mol, 106 pm
print('H2+ D0 kJ/mol', 2.65*96.485)
assert round(2.65*96.485) == 256
print('OK')
