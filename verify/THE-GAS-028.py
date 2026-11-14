# THE-GAS-028 検証：低圧での Z の傾きから van der Waals 定数 a, b を逆算，ボイル温度，Z=1 に戻る圧力
from sympy import Rational as Q, symbols, solve, N, nsolve
R = Q(831,100)*1000                      # Pa L/(K mol)
a0, b0 = Q(230,1)*1000, Q(430,10000)     # Pa L^2/mol^2, L/mol
def slope(T): return (b0 - a0/(R*T))/(R*T)        # P→0 での dZ/dP（vdW で厳密）
print('slope', N(slope(200), 6), N(slope(400), 6))
sg200, sg400 = Q(-574,100)*Q(1,10**8), Q(-788,100)*Q(1,10**9)   # 問題で与える値
assert abs(slope(200)/sg200 - 1) < 1e-3 and abs(slope(400)/sg400 - 1) < 1e-3
# 数値微分で vdW の Z(P) の P→0 の傾きを確認
V = symbols('V', positive=True)
for T in (200, 400):
    Pp = 10
    Vs = nsolve(R*T/(V - b0) - a0/V**2 - Pp, V, R*T/Pp)
    print(T, 'numerical slope', N((Pp*Vs/(R*T) - 1)/Pp, 5))
# 学生の手順：slope*RT = b - a/(RT)
a, b = symbols('a b')
s = solve([sg200*R*200 - (b - a/(R*200)), sg400*R*400 - (b - a/(R*400))], [a, b])
print('a', N(s[a]), 'b', N(s[b]))
assert abs(s[a] - 2.3e5) < 0.02e5 and abs(s[b] - 0.043) < 0.001
TB = s[a]/(R*s[b]); print('TB', N(TB)); assert abs(TB - 640) < 10
Tc = Q(8,27)*TB; print('Tc', N(Tc))
# 300 K, 800 K の低圧での符号
for T in (300, 800): print(T, 'slope sign', N(b0 - a0/(R*T)))
# 200 K で Z が再び 1 になるモル体積と圧力（a0, b0 で）
V1 = a0*b0/(a0 - R*200*b0); P1 = R*200/V1
print('V1', N(V1), 'P1', N(P1)); assert abs(P1 - 2.7e7) < 0.05e7
# 学生値 a,b でも同じ有効数字になるか
V1s = s[a]*s[b]/(s[a] - R*200*s[b]); print('P1(student)', N(R*200/V1s))
# vdW 式で Z=1 の点を直接確認
print('check P at V1', N(R*200/(V1 - b0) - a0/V1**2), 'RT/V1', N(R*200/V1))
print('OK')
