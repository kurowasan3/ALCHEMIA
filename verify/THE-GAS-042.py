# THE-GAS-042 検証：液体を使った気圧計（液柱の高さと蒸気圧）
from sympy import Rational as Q, N
P = 760; rHg = Q(136,10)                       # mmHg, g/cm3
h = lambda pv, rho: (P - pv)*rHg/rho          # mm
# (1) 水
h0 = h(0, 1); h1 = h(Q(238,10), 1)
print('(1) 無視', N(h0,5), ' 考慮', N(h1,5), ' 差', N(h0-h1,4), 'mm')
assert round(float(h0), -1) == 10340 and round(float(h1), -1) == 10010 and round(float(h0-h1)) == 324
# (2) エタノール（密度 0.79, 蒸気圧 59.0 mmHg）
he = h(59, Q(79,100)); print('(2) エタノール', N(he,5), 'mm; 蒸気圧の効果', N(59*rHg/Q(79,100),4))
assert he > h1 and round(float(he), -1) == 12070
# (3) 1.00 mol/kg NaCl 水溶液（密度 1.04）：ラウールの法則（粒子 2 mol）
xw = Q(1000,18)/(Q(1000,18)+2); pvs = Q(238,10)*xw
hs = h(pvs, Q(104,100))
print('(3) x_w', N(xw,5), ' p =', N(pvs,4), ' 降下', N(Q(238,10)-pvs,3), ' h =', N(hs,5))
print('    蒸気圧効果 +', N((Q(238,10)-pvs)*rHg,3), 'mm  密度効果', N(h(pvs,1)*(1/Q(104,100)-1),4), 'mm')
assert hs < h1 and round(float(hs), -1) == 9640
# (4) 密度 0.80 の液体で 11.5 m → 蒸気圧
pv = P - Q(11500)*Q(80,100)/rHg; print('(4) p =', N(pv,4), 'mmHg')
assert round(float(pv)) == 84
print('OK')
