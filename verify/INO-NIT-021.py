# INO-NIT-021 検証：肥料の成分表示（N・P2O5・K2O の質量%）と配合
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from sympy import Rational as Q, N
from chemlib import formula_mass, check_equation
m = dict(H=Q(1), C=Q(12), N=Q(14), O=Q(16), P=Q(31), S=Q(32), Cl=Q(355,10), K=Q(39), Ca=Q(40))
M = lambda f: formula_mass(f, m)
AS = M('(NH4)2SO4'); SP = M('Ca(H2PO4)2') + 2*M('CaSO4'); KCl = M('KCl'); UR = M('CO(NH2)2'); TSP = M('Ca(H2PO4)2')
P2O5, K2O = M('P2O5'), M('K2O')
print('式量', AS, SP, KCl, UR, TSP, P2O5, K2O)
fN_AS = 2*m['N']/AS
fP_SP = P2O5/SP           # Ca(H2PO4)2 1 つに P 2 個 = P2O5 1 個分
fK = (K2O/2)/KCl
fN_UR = 2*m['N']/UR
fP_TSP = P2O5/TSP
for n, v in [('N(硫安)', fN_AS), ('P2O5(過石)', fP_SP), ('K2O(KCl)', fK), ('N(尿素)', fN_UR), ('P2O5(重過石)', fP_TSP)]:
    print(n, N(100*v, 4), '%')
# (2) 8-8-8 を 100 g
a, b, c = 8/fN_AS, 8/fP_SP, 8/fK
print('(2) 硫安', N(a,4), '過石', N(b,4), 'KCl', N(c,4), '計', N(a+b+c,4), '増量材', N(100-a-b-c,4))
# (3) 15-15-15
a2, b2, c2 = 15/fN_AS, 15/fP_SP, 15/fK
print('(3) 15-15-15 硫安系 合計', N(a2+b2+c2, 4))
a3, b3 = 15/fN_UR, 15/fP_TSP
print('    尿素・重過石・KCl 合計', N(a3+b3+c2, 4), N(a3,4), N(b3,4), N(c2,4))
# 最大の等量表示 x-x-x（硫安系）
xmax = 100/(1/fN_AS + 1/fP_SP + 1/fK)
print('    硫安系の最大 x', N(xmax, 4))
assert a2+b2+c2 > 100 and a3+b3+c2 < 100
for e in ['Ca3(PO4)2 + 2H2SO4 -> Ca(H2PO4)2 + 2CaSO4', 'Ca3(PO4)2 + 4H3PO4 -> 3Ca(H2PO4)2']:
    assert check_equation(e) is None, e
print('OK')
