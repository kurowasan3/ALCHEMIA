# THE-EQL-150 verification: catalyst vs intermediate, V2O5 cycle, turnover number
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from collections import Counter
for e in ['H2O2 + I- -> IO- + H2O','H2O2 + IO- -> I- + H2O + O2','2H2O2 -> 2H2O + O2','SO2 + V2O5 -> SO3 + 2VO2','4VO2 + O2 -> 2V2O5','2SO2 + O2 -> 2SO3','MnO2 + 4HCl -> MnCl2 + Cl2 + 2H2O','CH3COOC2H5 + NaOH -> CH3COONa + C2H5OH','CH3COOC2H5 + H2O <=> CH3COOH + C2H5OH']:
    assert check_equation(e) is None, e
# sum of steps: 2*(SO2+V2O5->SO3+2VO2) + (4VO2+O2->2V2O5)
L=Counter({'SO2':2,'V2O5':2,'VO2':4,'O2':1}); R=Counter({'SO3':2,'VO2':4,'V2O5':2})
net_L=L-R; net_R=R-L; print(net_L,net_R); assert net_L==Counter({'SO2':2,'O2':1}) and net_R==Counter({'SO3':2})
nV=1820/(2*50.9+5*16.0); nS=6.40e6/64.0; print(nV,nS/nV); assert round(nS/nV,-3)==10000
print('OK')
