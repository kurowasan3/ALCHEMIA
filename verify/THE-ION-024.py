# THE-ION-024 検証：密閉容器中の CO2 と NaOH 水溶液（ヘンリーの法則＋炭酸の電離平衡＋物質収支）
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, symbols, nsolve, N, log
for e in ['CO2 + 2NaOH -> Na2CO3 + H2O', 'Na2CO3 + CO2 + H2O -> 2NaHCO3', 'CO2 + NaOH -> NaHCO3',
          'CO2 + 2OH- -> CO3^2- + H2O', 'CO3^2- + CO2 + H2O -> 2HCO3-', 'CO2 + H2O -> H+ + HCO3-',
          'HCO3- -> H+ + CO3^2-', 'CO3^2- + H2O -> HCO3- + OH-']:
    r = check_equation(e); assert r is None, (e, r)
R = Q(83,10)*1000; T = 300; RT = R*T            # 2.49e6
kH = Q(33,10)*Q(1,10**7)                        # mol/(L Pa)   (1.0e5 Pa で 3.3e-2 mol/L)
K1, K2, Kw = Q(45,10)*Q(1,10**7), Q(47,10)*Q(1,10**11), Q(1,10**14)
Vg, Vl, nNa = 1, 1, Q(1,100)
p, h = symbols('p h', positive=True)
def full(nC, guess):
    co2 = kH*p; hco3 = K1*co2/h; co3 = K2*hco3/h; oh = Kw/h
    e1 = p*Vg/RT + (co2 + hco3 + co3)*Vl - nC                 # 炭素の収支
    e2 = nNa/Vl + h - hco3 - 2*co3 - oh                        # 電荷の収支
    s = nsolve([e1, e2], [p, h], guess, prec=30)
    return s[0], s[1]
# (2) CO2 1.50e-2 mol
coef = 1/RT*Vg + kH*Vl; print('coef', N(coef))
pA_approx = Q(5,1000)/coef; hA_approx = K1*kH*pA_approx/Q(1,100)
print('approx p', N(pA_approx), 'pH', N(-log(hA_approx,10)))
pA, hA = full(Q(15,1000), [pA_approx, hA_approx])
print('exact p', N(pA), 'pH', N(-log(hA,10)))
assert abs(pA - 6.8e3) < 0.05e3 and abs(N(-log(hA,10)) - 7.0) < 0.05
# (3) CO3^2- の割合
print('[CO3]/[HCO3] =', N(K2/hA))
assert K2/hA < 1e-3
# (4) CO2 4.0e-3 mol（NaOH の半分未満）
oh0 = Q(2,1000); co3_0 = Q(4,1000)
pB, hB = full(Q(4,1000), [Q(1,100), Kw/oh0])
print('case B p', N(pB), 'pH', N(-log(hB,10)))
assert pB < Q(1,10)
# 近似の評価（加水分解を考慮した粗い見積り）
hco3_est = (Kw/K2)*co3_0/oh0; co2_est = (Kw/oh0)*hco3_est/K1
print('estimate p', N(co2_est/kH))
# (5) 逆算：p = 3.0e3 Pa のとき注入した CO2
nC5 = Q(1,100) + coef*3000; print('n injected', N(nC5))
p5, h5 = full(nC5, [3000, 2e-7]); print('check exact p', N(p5))
assert abs(nC5 - 1.2e-2) < 0.05e-2 and abs(p5 - 3000) < 30
# (5') 注入量が 7.0e-3 mol（NaOH の半分〜全量の間）のとき：緩衝域で p はごく小さい
pC, hC = full(Q(7,1000), [1, 1e-10]); print('7.0e-3 mol: p', N(pC), 'pH', N(-log(hC,10)))
print('OK')
