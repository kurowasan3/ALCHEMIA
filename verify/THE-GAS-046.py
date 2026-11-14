# THE-GAS-046 検証：容積・温度のわからない反応管をヘリウムで校正してアンモニアの生成量を求める
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, Eq, solve, Rational as Q
R, T, VA = 8.3e3, 300, 2.49
a = VA/(R*T)                       # A の「圧力あたりの物質量」mol/Pa
assert abs(a - 1.0e-6) < 1e-12
# 真の B（例）：VB, TB を適当に与えても結論が同じになることを確認
ex = 0
for VB, TB in ((2.988, 600), (4.98, 1000), (3.735, 750)):
    b = VB/(R*TB)
    P1 = 2.00e5; nHe = a*P1
    P2 = nHe/(a+b)
    if abs(P2-1.25e5) > 1: continue
    ex += 1
    print('B の例 VB=%.3f L, TB=%d K → P2 = %.0f Pa' % (VB, TB, P2))
assert ex == 3
# 問題文の値から
P1, P2 = 2.00e5, 1.25e5
nHe = a*P1; nA = a*P2; nB = nHe - nA
k = nHe/P2
print('(1) He 全体', nHe, ' B に入った', nB, ' k =', k)
assert round(nB, 4) == 0.075
P3, P4 = 2.40e5, 1.30e5
n0 = a*P3; nN2, nH2 = n0/4, 3*n0/4
n_after = k*P4
x = n0 - n_after                     # N2 + 3H2 -> 2NH3 で全物質量は NH3 の生成量だけ減る
print('(2) n0 =', n0, ' 反応後 =', n_after, ' NH3 =', x)
assert abs(x - 0.032) < 1e-9
pNH3 = P4*x/n_after; conv = (x/2)/nN2
print('(3) p(NH3) =', pNH3, ' N2 の反応率 =', conv)
assert round(pNH3/1e4, 1) == 2.0 and round(conv*100) == 27
# (4) B を高温にすると真の k は小さい → 古い k を使うと n_after を過大，NH3 を過小に見積もる
b_new = 0.075/1.25e5 * 0.5           # 例：TB が 2 倍
k_true = a + b_new
x_true_for_same_P = n0 - k_true*P4
x_est = n0 - k*P4
print('(4) 真の NH3', x_true_for_same_P, ' 古い校正での見積もり', x_est)
assert x_est < x_true_for_same_P
assert check_equation('N2 + 3H2 <=> 2NH3') is None
print('OK')
