# THE-RDX-030 検証：アスコルビン酸とアスコルビン酸ナトリウムの混合物（中和滴定とヨウ素滴定の組み合わせ）
import sys, math
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from sympy import symbols, nsolve, solve, log
from chemlib import check_equation
K1, K2 = 8.0e-5, 2.0e-12
cN, cI = 0.0500, 0.0500
V1, V2 = 8.00e-3, 12.0e-3
a = cN*V1                 # H2A（第1段だけ中和）
b = cI*V2 - a             # NaHA
print('a,b', a, b); assert abs(a-4.00e-4) < 1e-12 and abs(b-2.00e-4) < 1e-12
H, C, O, Na = 1.0, 12, 16, 23
M_H2A = 6*C+8*H+6*O; M_NaHA = 6*C+7*H+6*O+Na
print('M', M_H2A, M_NaHA, 'g per 100 mL', a*10*M_H2A, b*10*M_NaHA)
# (1) フェノールフタレイン変色域での第2段の進み
for pH in (8.0, 9.0, 9.8):
    print('pH', pH, '[A2-]/[HA-]=', K2/10**-pH)
# 第1当量点（NaHA のみ）の pH ≒ (pK1+pK2)/2
print('pH eq1 approx', (-math.log10(K1) - math.log10(K2))/2)
# 第1当量点の厳密計算（c = 6.00e-4 mol / (10+V) mL）: 電荷収支
Vtot = 10.0e-3 + (a + 0)/cN   # 10 mL + 8 mL
c = (a+b)/Vtot
x = symbols('x', positive=True)
# [H+] from charge balance with total c, Na = c (NaHA 全量)
Kw = 1e-14
def f(h): return h + c - c*(K1*h + 2*K1*K2)/(h*h + K1*h + K1*K2) - Kw/h
lo, hi = 1e-12, 1e-4
for _ in range(200):
    mid = math.sqrt(lo*hi)
    if f(mid) > 0: hi = mid
    else: lo = mid
print('exact pH eq1', -math.log10(mid))
assert 7.7 < -math.log10(mid) < 8.0
# (3) V3 = (2a+b)/cN
V3 = (2*a + b)/cN*1000
print('V3 mL', V3); assert abs(V3-20.0) < 1e-9
# (4) もとの溶液（10 mL 中）の pH
CA, CB = a/10e-3, b/10e-3
y = symbols('y', positive=True)
sol = nsolve((CB + y)*y/(CA - y) - K1, y, 1e-4)
print('buffer exact pH', -math.log10(float(sol)), 'approx', -math.log10(K1) + math.log10(CB/CA))
assert abs(-math.log10(float(sol)) - 3.80) < 0.01
# (5) メチルオレンジ終点 pH 4.4：[HA-]/[H2A] = K1/[H+]
Hp = 10**-4.4
ratio = K1/Hp
print('[H+] at 4.4', Hp, 'ratio', ratio)
ratio_r = K1/4.0e-5     # 与えた log2 で 10^-4.4 = 4.0e-5
tot = a + b
HA_end = tot*ratio_r/(1+ratio_r)
nOH = HA_end - b
print('NaOH mol', nOH, 'V mL', nOH/cN*1000)
assert abs(nOH/cN*1000 - 4.00) < 1e-9
# H+ の消費の寄与（溶液中の自由な H+ の変化）は無視できるか
print('free H+ change (mol, ~18mL)', (10**-3.8 - 4e-5)*0.014)
for eq in ['C6H8O6 + I2 -> C6H6O6 + 2H+ + 2I-', 'C6H7O6- + I2 -> C6H6O6 + H+ + 2I-',
           'C6H8O6 + OH- -> C6H7O6- + H2O', 'H+ + OH- -> H2O', 'I2 + 2S2O3^2- -> 2I- + S4O6^2-']:
    rr = check_equation(eq); print(eq, rr); assert rr is None
print('OK')
