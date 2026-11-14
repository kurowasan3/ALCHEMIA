# THE-THC-064 検証：ダイヤモンド・黒鉛の原子化エネルギーと平均結合エネルギー
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import sqrt
print(check_equation('C + O2 -> CO2'))
CO=804; OO=498
for name,dH,nb in [('diamond',395.4,2),('graphite',393.5,1.5)]:
    at=2*CO-OO-dH
    print(name,'atomization',at,'bond',at/nb)
print('diamond d=',sqrt(3)*0.357/4,' graphite d=',0.246/sqrt(3))
print(2*357.3, 1.5*477.7)
