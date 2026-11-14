# INO-TRN-043 シアン化法による金の回収 検証
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['4Au + 8CN- + O2 + 2H2O -> 4[Au(CN)2]- + 4OH-','2[Au(CN)2]- + Zn -> 2Au + [Zn(CN)4]^2-','CN- + H+ -> HCN','NaCN + HCl -> NaCl + HCN']:
    r=check_equation(e); print(e,r); assert r is None
nAu=5.0/197
NaCN=2*nAu*49; O2=nAu/4; Zn=nAu/2*65.4
print(nAu,NaCN,O2,Zn)
assert round(NaCN,1)==2.5 and round(O2*1e3,1)==6.3 and round(Zn,2)==0.83
# 電子の収支: Au 4 mol が 4 e- を O2 1 mol (4 e-) に渡す
print('OK')
