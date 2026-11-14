# THE-BND-058 検証：ケイ素単結晶球によるアボガドロ定数の決定（X線結晶密度法）
from sympy import Rational as R, nsimplify, N, symbols, solve, Eq
m28, m29, m30 = 27.976927, 28.976495, 29.973770
# 天然・濃縮の存在比（物質量の割合）
nat = (0.92223, 0.04685, 0.03092)
x29, x30 = 4.0e-5, 1.0e-6
enr = (1 - x29 - x30, x29, x30)
M = lambda x: x[0]*m28 + x[1]*m29 + x[2]*m30
Mn, Me = M(nat), M(enr)
print('(2) M 濃縮 =', Me, '-> 27.9770', ' 天然 =', Mn)
assert abs(Me - 27.97697) < 1e-5
# (1)(3) N_A = 8 M V / (m a^3)
NA, Mm, V, m, a = symbols('N_A M V m a', positive=True)
# 球中の原子数 = m N_A / M = 単位格子数 V/a^3 × 8
sol = solve(Eq(m*NA/Mm, 8*V/a**3), NA)[0]
print('(1) N_A =', sol)
V_, m_, a3 = 431.02, 1000.00, 0.16019e-21   # cm^3, g, cm^3（0.16019 nm^3）
NA_val = 8*27.9770*V_/(m_*a3)
print('(3) N_A =', NA_val)
assert abs(NA_val/1e23 - 6.022) < 0.0006
# 密度の確認
print('   密度 =', m_/V_, 'g/cm3（Si の 2.33 に近い）')
# (4) 少量同位体の存在比に相対 1% のずれ（どちらも同じ向き，28Si が補う）
dMn = 0.01*(nat[1]*(m29-m28) + nat[2]*(m30-m28))
dMe = 0.01*(enr[1]*(m29-m28) + enr[2]*(m30-m28))
print('(4) 天然 dM =', dMn, ' 相対', dMn/Mn, '  濃縮 dM =', dMe, ' 相対', dMe/Me)
assert abs(dMn/Mn - 3.9e-5) < 0.1e-5 and abs(dMe/Me - 1.5e-8) < 0.1e-8
# (5) 空孔の割合 f：実際の原子数は 8(1-f) 個/格子．8 個として計算した N_A は 1/(1-f) 倍
f = 3e-8
print('(5) 計算値/真の値 =', 1/(1-f), ' 相対ずれ', 1/(1-f)-1, '> 目標 2e-8')
assert 1/(1-f) - 1 > 2e-8
print('OK')
