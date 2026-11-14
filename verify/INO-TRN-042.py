# INO-TRN-042 金・銀・銅合金の分析 検証
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['Ag + 2HNO3 -> AgNO3 + NO2 + H2O','Cu + 4HNO3 -> Cu(NO3)2 + 2NO2 + 2H2O','AgNO3 + HCl -> AgCl + HNO3',
          'Au + HNO3 + 4HCl -> H[AuCl4] + NO + 2H2O','Cu^2+ + 4NH3 -> [Cu(NH3)4]^2+']:
    r=check_equation(e); print(e,r); assert r is None
m=1.000; Au=0.750; AgCl=0.166
Ag=AgCl*108/143.5; Cu=m-Au-Ag
print(Ag,Cu); assert round(Ag*100,1)==12.5 and round(Cu*100,1)==12.5
print('karat',24*Au/m)
nNO2=Ag/108+2*Cu/63.5; print('NO2',nNO2); assert round(nNO2*1e3,2)==5.10
print('OK')
