# INO-NIT-027 溶媒としての液体アンモニア 検証
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['2NH3 -> NH4+ + NH2-','2H2O -> H3O+ + OH-','2Na + 2NH3 -> 2NaNH2 + H2','2Na + 2H2O -> 2NaOH + H2',
          'NH4Cl + NaNH2 -> NaCl + 2NH3','HCl + NaOH -> NaCl + H2O','CH3COOH + NH3 -> CH3COO- + NH4+',
          '2NaN3 -> 2Na + 3N2','10Na + 2KNO3 -> K2O + 5Na2O + N2']:
    r=check_equation(e); print(e,r); assert r is None
nN2=1.5+0.1; V=nN2*8.31e3*300/1.01e5; print(V); assert round(V,1)==39.5
print('KNO3 g',0.2*101)
print('OK')
