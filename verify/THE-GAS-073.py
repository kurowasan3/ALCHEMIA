# THE-GAS-073 気体分子の表面への衝突頻度と単分子層ができる時間
from math import sqrt, pi, log10
R = 8.31; NA = 6.02e23; T = 300; pi_ = 3.14
M = 28.0e-3
v = sqrt(8*R*T/(pi_*M)); print('<v> N2 =', v)
P = 1.0e5
rho = P*NA/(R*T); print('rho =', rho)
Z = rho*v/4; print('Z =', Z)
sites = 1.0e19
t1 = sites/Z; print('t (1e5 Pa) =', t1)
# 1 時間以上清浄に保つ圧力
P_need = P*t1/3600; print('P for 1 h =', P_need)
# 1.0e-7 Pa
print('t at 1e-7 Pa =', t1*P/1e-7, 's =', t1*P/1e-7/3600, 'h')
# H2 との比較
print('Z(H2)/Z(N2) =', sqrt(28.0/2.0))
