# THE-ION-110 水酸化亜鉛の沈殿が存在する pH の範囲と濃度依存
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import log10, sqrt
Ksp=2.0e-17; beta=2.5e17; Kw=1e-14
K4=Ksp*beta
print('K4 =',K4)
for c in [1e-2,1e-1,1e-3]:
    oh1=sqrt(Ksp/c); oh2=sqrt(c/K4)
    pH1=14+log10(oh1); pH2=14+log10(oh2)
    print(c,'start pH %.3f end pH %.3f width %.2f'%(pH1,pH2,pH2-pH1),
          'check: complex at start',K4*oh1**2,'Zn2+ at end',Ksp/oh2**2)
# 有効数字（log2=0.30 使用）
print('c=0.010: pOH start=(15-0.30)/2=',(15-0.30)/2,'pOH end=(2+0.30... )')
print(check_equation('Zn(OH)2 + 2OH- -> [Zn(OH)4]^2-'))
print(check_equation('Zn^2+ + 2OH- -> Zn(OH)2'))
