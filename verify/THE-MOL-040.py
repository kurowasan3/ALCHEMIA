# THE-MOL-040 検証：当量とデュロン・プティの法則による原子量の決定
from sympy import Rational as Q, N
O = Q(16)
# 金属A：1.00 g -> 酸化物 1.89 g
mA, oxA = Q('1.00'), Q('1.89')
EA = 8*mA/(oxA-mA)             # 酸素 8.00 g と化合する質量（=X/n）
print('E_A =', N(EA, 6))       # 8.988...
approxA = 26/Q('0.90')
nA = round(approxA/EA); XA = nA*EA
print('DP approx A =', N(approxA, 4), ' n =', nA, ' X =', N(XA, 5))
assert round(float(EA), 2) == 8.99 and nA == 3 and round(float(XA), 1) == 27.0
# 金属B：二種の酸化物
mB = Q('1.00'); P, Qx = Q('1.252'), Q('1.126')
EP = 8*mB/(P-mB); EQ = 8*mB/(Qx-mB)
approxB = 26/Q('0.385')
nP = round(approxB/EP); nQ = round(approxB/EQ)
print('E_P', N(EP, 5), 'E_Q', N(EQ, 5), 'DP', N(approxB, 4), nP, nQ, N(nP*EP, 5), N(nQ*EQ, 5))
assert (nP, nQ) == (2, 1) and round(float(nP*EP), 1) == 63.5 and round(float(nQ*EQ), 1) == 63.5
# 組成式：P は B:O = 1:1（BO），Q は 2:1（B2O）
# ドルトン流（すべて1:1）での見かけの原子量
print('Dalton A =', N(16*mA/(oxA-mA), 4), ' B(P) =', N(16*mB/(P-mB), 4), ' B(Q) =', N(16*mB/(Qx-mB), 4))
assert round(float(16*mA/(oxA-mA)), 1) == 18.0
assert round(float(16/(P-1)), 1) == 63.5 and round(float(16/(Qx-1)), 0) == 127
# 実在金属との整合（Al, Cu）
print('Al2O3/2Al =', 101.96/53.96, ' CuO/Cu =', 79.55/63.55, ' Cu2O/2Cu =', 143.09/127.09)
print('OK')
