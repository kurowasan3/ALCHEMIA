# THE-ELC-086 検証：仕切りのある電解槽での水の電気分解と各室の pH
import sys, math; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['2H2O -> O2 + 4H+ + 4e-','2H2O + 2e- -> H2 + 2OH-','2H2O -> 2H2 + O2','H+ + OH- -> H2O']:
    assert check_equation(e) is None, e
F=9.65e4; I=0.193; t=500
ne=I*t/F; print('e- mol',ne); assert abs(ne-1.00e-3)<1e-12
V=0.100
H=ne/V; OH=ne/V
pHa=-math.log10(H); pHc=14+math.log10(OH)
print('陽極室 pH',pHa,'陰極室 pH',pHc)
assert round(pHa,1)==2.0 and round(pHc,1)==12.0
print('O2',ne/4,'H2',ne/2,'分解した水',ne/2)
assert abs(ne/4-2.5e-4)<1e-12 and abs(ne/2-5.0e-4)<1e-12
# 混合：H+ と OH- が同物質量で中和 → 中性
assert abs(H*V-OH*V)<1e-15
# 陰極室で生じた OH- の半分が中和されたとき（陽極室の液を 50 mL だけ加えた場合）の pH
rem=ne-ne/2; pH_mix=14+math.log10(rem/(V+0.050)); print('50 mL 混合',pH_mix)
assert round(pH_mix,1)==11.5
print('OK')
