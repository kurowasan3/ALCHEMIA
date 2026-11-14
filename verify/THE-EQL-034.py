# THE-EQL-034 検証：ヘモグロビン（全か無かモデル）とミオグロビンの酸素飽和度
from sympy import Rational as Q, nsolve, symbols, N, log
import sys
P50_Hb, P50_Mb = Q(36,10), Q(4,10)          # kPa
Y_Hb = lambda p, P50=P50_Hb: p**4/(p**4+P50**4)
Y_Mb = lambda p: p/(p+P50_Mb)
pl, pr, pa = 12, Q(48,10), Q(24,10)
vals = {k: N(f(p), 6) for k, f, p in [
    ('Hb12', Y_Hb, pl), ('Hb4.8', Y_Hb, pr), ('Hb2.4', Y_Hb, pa),
    ('Mb12', Y_Mb, pl), ('Mb4.8', Y_Mb, pr), ('Mb2.4', Y_Mb, pa)]}
print(vals)
rel = {'Hb_rest': N(Y_Hb(pl)-Y_Hb(pr), 4), 'Hb_act': N(Y_Hb(pl)-Y_Hb(pa), 4),
       'Mb_rest': N(Y_Mb(pl)-Y_Mb(pr), 4), 'Mb_act': N(Y_Mb(pl)-Y_Mb(pa), 4)}
print(rel)
assert abs(vals['Hb12']-0.99) < 0.005 and abs(vals['Hb4.8']-0.76) < 0.005 and abs(vals['Hb2.4']-0.165) < 0.002
assert abs(rel['Hb_act']-0.83) < 0.005 and abs(rel['Mb_act']-0.11) < 0.005 and abs(rel['Hb_rest']-0.23) < 0.005
# (3) 血液 100 mL：Hb 14 g，M = 6.4e4，1分子に O2 4個．活動時の放出量
nHb = Q(14)/64000; nsite = 4*nHb
nrel = nsite*(Y_Hb(pl)-Y_Hb(pa))
V = nrel*Q(83,10)*1000*310/10**5*1000   # mL（R=8.3e3 Pa L/(K mol), 37℃, 1.0e5 Pa）
print('nrel', N(nrel), 'mol  V', N(V), 'mL')
assert abs(N(nrel)-7.2e-4) < 0.06e-4 and abs(V-19) < 0.5
# (4) 積み込んだ O2 の 60 % を放出する組織の酸素分圧
p = symbols('p', positive=True)
sol = nsolve(Y_Hb(p) - Q(4,10)*Y_Hb(pl), p, 3)
print('p(60%放出)', sol); assert abs(sol-3.24) < 0.01
# (5) ボーア効果：P50 = 4.4 kPa になったとき活動筋での放出割合
r2 = N(Y_Hb(pl) - Y_Hb(pa, Q(44,10)), 4)
y2 = N(Y_Hb(pa, Q(44,10)), 4)
print('bohr', y2, r2); assert abs(r2-0.91) < 0.005
# 肺での積み込みは P50 が変わっても（肺では元の P50 で積み込む）0.99
print('OK')
