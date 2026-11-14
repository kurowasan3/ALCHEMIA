# THE-RDX-053 verification: Cr by BaCrO4 gravimetry, Mn by Fe2+ titration minus Cr share
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import log10,floor
def sig(v,n): return round(v,-int(floor(log10(abs(v))))+n-1)
for e in ['Cr2O7^2- + 14H+ + 6e- -> 2Cr^3+ + 7H2O','MnO4- + 8H+ + 5e- -> Mn^2+ + 4H2O','Fe^2+ -> Fe^3+ + e-','Cr2O7^2- + 2Ba^2+ + H2O -> 2BaCrO4 + 2H+','2CrO4^2- + 2H+ <=> Cr2O7^2- + H2O','CH3COO- + H+ -> CH3COOH','2MnO4- + 10Cl- + 16H+ -> 2Mn^2+ + 5Cl2 + 8H2O']:
    assert check_equation(e) is None, e
Mb=137.3+52.0+4*16.0
# forward design check
cr=0.170*2.50/52.0; mn=0.0200*2.50/54.9
print('BaCrO4',cr/5*Mb,'Fe mL',((cr*20/250)*3+(mn*20/250)*5)/0.100*1000)
nb=0.414/Mb; crp=nb*5*52.0/2.50*100
fe=0.100*23.3e-3; mn20=(fe-3*nb*20/50)/5; mnp=mn20*250/20*54.9/2.50*100
wrong=fe/5*250/20*54.9/2.50*100
print(crp,mnp,wrong); assert sig(crp,2)==17 and sig(mnp,2)==2.0 and sig(wrong,2)==13
print('OK')
