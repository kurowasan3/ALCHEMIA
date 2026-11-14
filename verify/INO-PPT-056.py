# INO-PPT-056 連続変化法（ジョブ法）：有限の安定度定数と極大の位置
from sympy import symbols, Rational as R, solve, nsolve, diff, simplify, Eq
import numpy as np
from scipy.optimize import brentq
C0 = 2.00e-3; eps = 1.50e3; L = 1.00
# (1) K 無限大：A(x)
def A_inf(x): return eps*L*min(x/2, 1-x)*C0
print('ideal peak at 2/3:', A_inf(2/3), ' A(0.5)=', A_inf(0.5), ' A(0.8)=', A_inf(0.8))
# (2) x=2/3 で実測 0.800 → K
y = 0.800/(eps*L); Mt = C0/3; Lt = 2*C0/3
K = y/((Mt-y)*(Lt-2*y)**2); print('K =', K)
# (3) 有限 K で極大位置が 2/3 か数値確認
def yfun(x, K, C0):
    Mt, Lt = (1-x)*C0, x*C0
    f = lambda yy: yy - K*(Mt-yy)*(Lt-2*yy)**2
    return brentq(f, 0, min(Mt, Lt/2)*(1-1e-15))
for KK, CC in ((K, C0), (K, C0/10), (1e5, 1e-2)):
    xs = np.linspace(0.3, 0.95, 130001)
    ys = [yfun(x, KK, CC) for x in xs[::100]]
    i = int(np.argmax(ys)); xc = xs[::100][i]
    xs2 = np.linspace(xc-0.002, xc+0.002, 4001); ys2 = [yfun(x, KK, CC) for x in xs2]
    print(f'K={KK:.3g}, C0={CC}: x_max = {xs2[int(np.argmax(ys2))]:.5f}')
# 記号での証明：dy=0 のとき l=2m → Lt=2Mt
# (4) 薄めたとき：C0/10 での x=2/3 の A/A_ideal
for CC in (C0, C0/10, C0*10):
    yy = yfun(2/3, K, CC); print(f'C0={CC}: ratio y/(C0/3) = {yy/(CC/3):.3f}')
# ML3 の場合の極大位置 n/(n+1)
print('ML3 ideal peak', 3/4)
