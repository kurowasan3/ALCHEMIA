# THE-GAS-067 検証：圧力一定・不活性気体共存下での液体の蒸発
from sympy import *
R=8.31e3  # Pa L/(K mol)
P=94.0    # kPa
pv={30:10.4,40:17.9,50:29.5,60:47.0,70:72.3,80:108.3}
ne=0.100; nN=0.100
# (1) 40℃
pN=P-pv[40]; nv=nN*pv[40]/pN; liq=ne-nv
V=nN*R*313/(pN*1e3)
print('(1) nv=',nv,'liq=',liq,'frac liq',liq/ne,'V=',V)
# (2) 液体が消える：pv(T)=P*ne/(ne+nN)
print('(2) target pv=',P*ne/(ne+nN))
# (3) 80℃ のエタノール分圧（液体なし）と体積
V80=(ne+nN)*R*353/(P*1e3); print('(3) p_e=',P/2,'V80=',V80)
# (4) N2 0.200 mol
t=P*ne/(ne+0.200); print('(4) target pv=',t,' between 50(29.5) and 60(47.0)')
# 線形補間の目安
print(50+10*(t-29.5)/(47.0-29.5))
