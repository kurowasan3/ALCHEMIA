# THE-RDX-050 ヨウ素→アスコルビン酸→銀イオンの順に標定・定量 検証
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['C6H8O6 + I2 -> C6H6O6 + 2HI','C6H8O6 + 2Ag+ -> C6H6O6 + 2Ag + 2H+','I2 + I- -> I3-','C6H8O6 -> C6H6O6 + 2H+ + 2e-','Ag+ + e- -> Ag']:
    r=check_equation(e); print(e,r); assert r is None
cI2=0.508/254/0.200; print('I2',cI2); assert abs(cI2-1.00e-2)<1e-9
nasc=cI2*20.0e-3; cB=nasc/25.0e-3; print('B',cB); assert abs(cB-8.00e-3)<1e-9
nAg10=2*cB*10.5e-3; nAg=nAg10*10; mass=nAg*170; pct=mass/0.300*100
print(nAg,mass,pct); assert round(pct,1)==95.2
# B が空気酸化で 5 % 減っていた場合（実際 7.60e-3）：必要な体積は 10.5/0.95，計算値は過大
nAg_true=1.68e-4; V=nAg_true/2/(0.95*cB)   # 実際に要する体積 (L)
calc=2*cB*V*10*170/0.300*100
print('needed mL',V*1e3,'computed %',calc); assert calc>pct
print('OK')
