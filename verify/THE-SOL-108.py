# THE-SOL-108 検証：逆滴定（空試験）による有機層濃度と会合数の判定
from sympy import Rational as R, nsimplify
import math
# (1) 水層の分析濃度
Cw1 = 0.0200*3.10/20.0
print('Cw1 =', Cw1)            # 3.10e-3
# (2) 逆滴定：有機層10.0mL→NaOH 25.0mL に移し，10.0mL を HCl 0.0500 で滴定。空試験 20.40, 試料 17.20
nHCl_diff = 0.0500*(20.40-17.20)/1000      # 分取 10.0 mL 中の HA
nHA = nHCl_diff*25.0/10.0                   # 有機層 10.0 mL 中
Co1 = nHA/(10.0/1000)
print('Co1 =', Co1)            # 4.00e-2
# (3) 表のデータ
data = [(Co1, Cw1), (9.37e-3, 1.50e-3), (1.50e-3, 6.00e-4)]
for n in (1,2,3):
    print('n=',n, [round(co/cw**n, 4 if n==1 else 0) for co,cw in data])
D = [co/cw for co,cw in data]; print('D', [round(d,3) for d in D])
K2 = sum(co/cw**2 for co,cw in data)/3
print('K(n=2) mean', K2)
# (4) Cw=2.00e-3 での予測
Cw = 2.00e-3; Co = 4.16e3*Cw**2
print('pred Co', Co, 'D', Co/Cw)
# 表の値の作り方の確認：K=4.162e3
K = Co1/Cw1**2
print('K from run1', K, K*1.50e-3**2, K*6.00e-4**2)
