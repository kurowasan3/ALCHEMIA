# THE-SOL-105 検証：浸透圧から分子量・断片の数・D の分子量
from sympy import Rational, nsimplify
R = 8.31e3; T = 300.0
w = 1.50; V = 0.100
Pi0 = 1.50e3
MB = w*R*T/(Pi0*V); print('M_B =', MB)            # 2.49e4 -> 2.5e4
PiX = 6.0e3
ratio = PiX/Pi0; n = round(ratio) - 1; print('ratio', ratio, 'n =', n)
assert abs(ratio-4) < 1e-9 and n == 3
ME = 600
MD = MB - n*ME; print('M_D =', MD)                 # 2.31e4 -> 2.3e4
# 膜Y：平衡では E の濃度が内外で等しく浸透圧に寄与しない。D の数＝B の数
PiY_final = Pi0; print('Pi_Y final =', PiY_final)
# 分解直後（E がまだ膜を通っていない）の浸透圧
print('Pi immediately =', (1+n)*Pi0)
assert f'{MB:.2g}' == '2.5e+04' and f'{MD:.2g}' == '2.3e+04'
print('OK')
