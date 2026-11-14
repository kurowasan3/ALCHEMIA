# THE-MOL-052 検証：水素化ホウ素の分子式・電子の数・加水分解
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R
B, H = R('10.8'), R('1.00')
for x, y in [(4, 10), (6, 10)]:
    M = x*B + y*H
    print(f'B{x}H{y}', float(M), 'B%', float(x*B/M*100))
# 逆に：分子量と B% から x, y を求める（整数になるか）
for M, p in [(R('53.2'), R('81.2')), (R('74.8'), R('86.6'))]:
    x = M*p/100/B; y = (M - M*p/100)/H
    print('M', float(M), 'x', float(x), 'y', float(y))
# 電子の数：B2H6，B4H10
for name, nB, nH, term, BB, bridge in [('B2H6', 2, 6, 4, 0, 2), ('B4H10', 4, 10, 6, 1, 4)]:
    ve = 3*nB + nH
    used = 2*term + 2*BB + 2*bridge
    allsingle = 2*(term + BB + 2*bridge)
    print(name, 'valence e', ve, 'used', used, 'if all 2c', allsingle, 'short', allsingle - ve)
    assert ve == used
for e in ['B2H6 + 6H2O -> 2H3BO3 + 6H2', 'B4H10 + 12H2O -> 4H3BO3 + 11H2', 'B6H10 + 18H2O -> 6H3BO3 + 14H2', 'B2H6 + 2H- -> 2BH4-', 'B2H6 + 2NaH -> 2NaBH4']:
    r = check_equation(e); print(e, r); assert r is None
n = R('2.66')/(4*B + 10*H)
V = n*11*R('22.4')
print('H2 L', float(V))
assert f'{float(V):.3g}' == '12.3'
