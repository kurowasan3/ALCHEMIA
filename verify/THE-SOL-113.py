# THE-SOL-113 検証：溶解結晶法（90 ℃ と 10 ℃ の NaCl-KCl-水 溶解度図）の物質収支
from sympy import Rational as R, symbols, solve, nsimplify
N, K = symbols('N K')
# 90 ℃：A(0,53.0) B(28.0,34.0) C(38.5,0)  10 ℃：A'(0,31.0) B'(30.0,12.0) C'(35.8,0)
def line(p, q):
    (x1, y1), (x2, y2) = p, q
    return lambda n: y1 + (y2 - y1)*(n - x1)/(x2 - x1)
A, B, C = (R(0), R(53)), (R(28), R(34)), (R('38.5'), R(0))
A1, B1, C1 = (R(0), R(31)), (R(30), R(12)), (R('35.8'), R(0))
Kc = line(A1, B1)       # 10 ℃ KCl 飽和線 K(N)
# (1) B を 10 ℃ に急冷：N=28 <= 30 なので KCl 側
KT = Kc(R(28)); dep_K = 34 - KT
print('(1) T =', (28, float(KT)), 'KCl deposited', float(dep_K))
assert f'{float(dep_K):.3g}' == '20.7'
# NaCl が析出しないことの確認：T は N<=30 で A'B' 上
assert 28 < 30
# (2) T を 90 ℃ で濃縮：原点からの半直線 K/N = KT/28 が BC に当たる点 U
s = KT/28
NU = solve(N - (C[0] + (B[0]-C[0])*(s*N - C[1])/(B[1]-C[1])), N)[0]
KU = s*NU
print('(2) U =', float(NU), float(KU))
W_U = 28*100/NU      # 水 100 g（最初）あたりの残りの水
print('   water at U', float(W_U), 'evaporated', float(100 - W_U))
W_B = KT*100/34      # B に達したときの水（KCl は全部溶液中）
NaCl_sol = 28*W_B/100
dep_N = 28 - NaCl_sol
print('   water at B', float(W_B), 'evaporated total', float(100 - W_B), 'NaCl deposited', float(dep_N))
assert f'{float(100 - W_U):.3g}' == '16.6'
assert f'{float(W_B):.1f}' == '39.0' and f'{float(dep_N):.3g}' == '17.1'
# (3) 定常運転：取り出した KCl, NaCl を原料で補う
frac = dep_K/(dep_K + dep_N)
print('(3) KCl fraction in ore', float(frac), 'water to add', float(100 - W_B))
assert f'{float(frac)*100:.3g}' == '54.8'
# (4) 低温の共飽和点が (26.0, 13.0) なら：B(28,34) は 10 ℃ で領域 III -> 水を加えて N<=26 に
w_add = 28*R(100)/26 - 100
print('(4) water to add', float(w_add))
assert f'{float(w_add):.2g}' == '7.7'
# 希釈後 (26, 34*100/107.69) は KCl 側で KCl のみ析出（K > 13）
print('   diluted point', 26, float(34*100/(100 + w_add)))
