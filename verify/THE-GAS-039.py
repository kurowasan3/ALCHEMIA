# THE-GAS-039 検証：-20 ℃ の氷 9.0 g に熱を加えたときの最終状態
from sympy import Rational as Q, N
m = Q(9); M = 18; n = Q(m, M)
c_ice, c_w = Q(21,10), Q(42,10)          # J/(g K)
Hf, Hv = 6000, 41000                     # J/mol
qA = m*c_ice*20; qB = n*Hf; qC = m*c_w*100; qD = n*Hv
print('段階 A,B,C,D [kJ]:', [N(x/1000,4) for x in (qA,qB,qC,qD)])
assert max(qA,qB,qC,qD) == qD
print('(1) D/B =', N(qD/qB,4), ' C/B =', N(qC/qB,4))
cum = [qA, qA+qB, qA+qB+qC, qA+qB+qC+qD]
print('累積 [kJ]:', [N(x/1000,5) for x in cum])
# (2) 5.0 kJ
Q2 = 5000
assert cum[1] < Q2 < cum[2]
t2 = (Q2-cum[1])/(m*c_w); print('(2) t =', N(t2,4), '℃')
assert round(float(t2)) == 43
# (3) 20.0 kJ
Q3 = 20000
assert cum[2] < Q3 < cum[3]
nv = (Q3-cum[2])/Hv; left = m - nv*M
print('(3) 蒸発', N(nv,4), 'mol', N(nv*M,4), 'g, 残り液体', N(left,4), 'g')
assert abs(float(left)-3.4) < 0.05
# (4) 沸点 80 ℃（蒸発熱は同じとする）
qC2 = m*c_w*80
nv4 = (Q3-qA-qB-qC2)/Hv; left4 = m - nv4*M
print('(4) 残り液体', N(left4,4), 'g（(3)より少ない）', '差', N(left-left4,3))
assert left4 < left and abs(float(left4)-3.0) < 0.05
print('OK')
