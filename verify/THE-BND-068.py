# THE-BND-068 検証：O2- の六方最密配列のすき間占有と密度
from sympy import Rational as Q, sqrt, N, nsimplify
# 六角柱（12頂点×1/6＋2面心×1/2＋内部3）＝6個
nO = 12*Q(1,6) + 2*Q(1,2) + 3
assert nO == 6
octa = nO; tetra = 2*nO
# Al2O3：Al = 2/3 × O
nAl = Q(2,3)*nO
print('Al per prism', nAl, 'occupancy', nAl/octa)
assert nAl == 4 and nAl/octa == Q(2,3)
# 六角柱の体積 = (3√3/2) a^2 c, c = √(8/3) a → 3√2 a^3
a = 1
V = Q(3,2)*sqrt(3)*a**2*sqrt(Q(8,3))*a
print('V/a^3 =', nsimplify(V))
assert nsimplify(V - 3*sqrt(2)) == 0
# 密度（√2 = 1.41 を用いる）
a_cm = 2.75e-8
Vcm = 3*1.41*a_cm**3
m = (4*27 + 6*16)/6.0e23
rho = m/Vcm
print('V =', Vcm, 'm =', m, 'rho =', rho)
assert round(rho, 1) == 3.9
# 厳密な √2 でも 2桁は同じ
print('rho(exact sqrt2) =', m/(3*2**0.5*a_cm**3))
assert round(m/(3*2**0.5*a_cm**3), 1) == 3.9
# ABO3：陽イオン2個で +6
pairs = {'Fe2+,Ti4+': 2+4, 'Li+,Nb5+': 1+5, 'Mg2+,Al3+': 2+3, 'Na+,Ti4+': 1+4}
ok = [k for k, v in pairs.items() if v == 6]
print(ok); assert ok == ['Fe2+,Ti4+', 'Li+,Nb5+']
# ZnO：Zn2+ 数 = O2- 数 = 四面体すき間の 1/2
assert Q(nO, tetra) == Q(1,2)
# 四面体すき間すべてに +1 → M2O
assert tetra/nO == 2
print('OK')
