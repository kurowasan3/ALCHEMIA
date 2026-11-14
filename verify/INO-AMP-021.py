# INO-AMP-021 verification: AlCl3 dimer, vapour density, hydrate decomposition
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
R=8.31e3
M1=6.1*R*523/1.0e5; M2=1.46*R*1100/1.0e5; print(M1,M2)
assert round(M1)==265 and round(M2)==133
Al,Cl,H,O=27.0,35.5,1.0,16.0
Mh=Al+3*Cl+6*(2*H+O); n=4.83/Mh; m=n/2*(2*Al+3*O); print(Mh,n,m); assert round(m,2)==1.02
for e in ['2Al + 3Cl2 -> 2AlCl3','2Al + 6HCl -> 2AlCl3 + 3H2','2[Al(H2O)6]Cl3 -> Al2O3 + 6HCl + 9H2O','AlCl3 + 3H2O -> Al(OH)3 + 3HCl','2AlCl3 -> Al2Cl6']:
    assert check_equation(e) is None, e
print('OK')
