# THE-SOL-064 検証：当量（質量変化）と価数（凝固点降下）を別々の情報から決める
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, N, symbols, solve
for e in ['Ba(NO3)2 + H2SO4 -> BaSO4 + 2HNO3']:
    r = check_equation(e); print(e, r or 'OK'); assert r is None
# 実在の Ba で数値を作った確認（Ba 137.33, N 14.007, O 15.999, S 32.06）
Ba, Nn, O, S = 137.33, 14.007, 15.999, 32.06
nit = Ba + 2*(Nn + 3*O); sul = Ba + S + 4*O
print('2.000 g Ba(NO3)2 -> BaSO4', 2.000*sul/nit)   # 1.786
# (1) 当量 E = X/n
w1, w2 = Q('2.000'), Q('1.786')
E = symbols('E')
Es = solve(w1/(E + 62) - w2/(E + 48), E)[0]
print('E =', N(Es, 6)); assert round(float(Es), 1) == 68.8
# (2) 完全電離での凝固点降下度（硝酸塩 2.000 g / 水 100 g, Kf = 1.85）
Kf = Q('1.85')
dT = {}
for n in (1, 2, 3):
    eq = w1/(Es + 62)                 # 正電荷（=負電荷）の物質量 = 当量の物質量
    particles = eq*(Q(1, n) + 1)       # 陽イオン eq/n 個 + NO3- eq 個
    dT[n] = Kf*particles/Q('0.100')
    print('n =', n, 'X =', N(n*Es, 5), 'dT =', N(dT[n], 5))
assert [round(float(dT[n]), 3) for n in (1, 2, 3)] == [0.566, 0.424, 0.377]
# 直接の式量計算との一致
for n in (1, 2, 3):
    X = n*Es; Mf = X + 62*n
    assert abs(N(Kf*w1/Mf*(1 + n)/Q('0.100') - dT[n])) < 1e-12
# (3) 測定値 0.420 K, 相対誤差 3% 以内
meas = Q('0.420')
ok = [n for n in (1, 2, 3) if abs(meas - dT[n])/dT[n] <= Q('0.03')]
print('consistent n:', ok, {n: N((meas - dT[n])/dT[n], 3) for n in dT}); assert ok == [2]
X = 2*Es; print('X =', N(X, 4)); assert round(float(X)) == 138
# 候補：Ga 69.7, Ba 137.3, Pb 207.2, Tl 204.4, Bi 209.0 → Ba
# 当量だけでは n=3 (X=206.5) で Pb/Tl に近く区別できない：E の不確かさ
for dw in (-0.0005, 0.0005):
    Ed = solve(w1/(E + 62) - (w2 + Q(dw))/(E + 48), E)[0]
    print('E range', N(Ed, 5), '3E =', N(3*Ed, 5))
# (4) 陰イオンの価数 a による識別力：ΔT ∝ (1/n + 1/a)
for a, name in ((1, 'NO3-'), (2, 'SO4^2-')):
    r = (Q(1, 2) + Q(1, a))/(Q(1, 3) + Q(1, a))
    print(name, 'ratio n=2/n=3 =', r, N(r, 4))
assert (Q(1,2)+1)/(Q(1,3)+1) == Q(9, 8) and (Q(1,2)+Q(1,2))/(Q(1,3)+Q(1,2)) == Q(6, 5)
print('OK')
