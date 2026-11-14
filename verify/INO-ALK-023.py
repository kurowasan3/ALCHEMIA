# INO-ALK-023 カルシウム化合物の識別 検証
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
eqs=['CaC2 + 2H2O -> Ca(OH)2 + C2H2','CaH2 + 2H2O -> Ca(OH)2 + 2H2','CaO + H2O -> Ca(OH)2',
     'Ca(ClO)2 + 4HCl -> CaCl2 + 2Cl2 + 2H2O','CaCO3 + 2HCl -> CaCl2 + H2O + CO2',
     'CaCO3 + H2SO4 -> CaSO4 + H2O + CO2','2C2H2 + 5O2 -> 4CO2 + 2H2O','2H2 + O2 -> 2H2O',
     'Ca(OH)2 + CO2 -> CaCO3 + H2O']
for e in eqs: r=check_equation(e); print(e, r); assert r is None
from fractions import Fraction as F
# 1 mol の A (CaC2) -> C2H2 1 mol -> O2 5/2 ; 1 mol の B (CaH2) -> H2 2 mol -> O2 1
ratio=F(5,2)/F(1); print(ratio); assert ratio==F(5,2)
# 同じ質量 4.2 g ずつなら: CaC2 64, CaH2 42
nA=F(42,10)/64; nB=F(42,10)/42
print(float(nA*F(5,2)), float(nB*1))
print('OK')
