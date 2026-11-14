# THE-ION-113 ホウ酸のルイス酸としての電離とポリオール（マンニトール）による見かけの酸性の増大
import sys, math; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
Ka=6.0e-10; beta=1.0e5; h=1e-9
Kp=lambda L: Ka*(1+beta*L**2)
print("Ka' at L=0.50: %.4g"%Kp(0.5))
f=lambda K: K/(K+h)
print('fraction anionic at pH9: none %.4f, L=0.50 %.6f'%(f(Ka),f(Kp(0.5))))
for K in [Ka,Kp(0.5)]:
    x=(-K+math.sqrt(K*K+4*K*0.10))/2; print('pH of 0.10 M %.3f'%(-math.log10(x)))
Kneed=0.999*h/0.001; Lmin=math.sqrt((Kneed/Ka-1)/beta)
print('Kneed %.4g Lfree min %.4f'%(Kneed,Lmin))
cB=0.050; comp=Ka*beta*Lmin**2/(h+Kp(Lmin)); Ltot=Lmin+2*cB*comp
print('complex frac %.4f  L total %.4f M  mol in 20.0 mL %.4e  g %.3f'%(comp,Ltot,Ltot*0.020,Ltot*0.020*182))
for e in ['B(OH)3 + H2O <=> [B(OH)4]^- + H+','B(OH)3 + OH- -> [B(OH)4]^-']:
    print(e,check_equation(e))
