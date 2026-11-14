# THE-ION-111 CO2 の水和を含めた見かけの電離定数と血液の pH
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import log10
Kh=10**-2.8; Ka=10**-3.6
Kapp=Ka*Kh/(1+Kh); print('Kapp',Kapp,'log',log10(Kapp),'approx Ka*Kh log',log10(Ka*Kh))
print('true H2CO3 vs acetic (10^-4.8): stronger?',Ka>10**-4.8,' apparent weaker?',Kapp<10**-4.8)
k=2.3e-7; P1=5.2e3; P2=7.8e3
c1=k*P1; h=10**-7.4; K37=10**-6.1
hco3=K37*c1/h; print('CO2 %.3g HCO3 %.3g ratio %.3g'%(c1,hco3,hco3/c1))
c2=k*P2; pH2=6.1+log10(hco3/c2); print('pH at P2 %.3f'%pH2, 'need HCO3 %.3g'%(hco3*1.5))
for e in ['CO2 + H2O <=> H2CO3','H2CO3 <=> H+ + HCO3-','CO2 + H2O <=> H+ + HCO3-']: print(e,check_equation(e))
