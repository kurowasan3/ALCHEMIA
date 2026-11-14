# THE-BND-099 verification: CaNi5 hexagonal prism, density, hydride composition, storage density
from math import log10,floor
def sig(v,n): return round(v,-int(floor(log10(abs(v))))+n-1)
Ca=12/6+2/2; Ni=12/2+6/2+6; assert (Ca,Ni)==(3,15)
NA=6.02e23; a=4.95e-8; c=3.94e-8
V=3*1.73/2*a**2*c; M=3*(40.1+5*58.7); rho=M/NA/V; print(V,rho); assert sig(rho,2)==6.6
nH=1.11/22.4*2*NA; perprism=nH*V; print(perprism); assert round(perprism)==15
mH=1.11/22.4*2.0; print(mH, mH/0.071); assert sig(mH/0.071,2)==1.4
gas=2.0e7/(8.31e3*298)*1e-3*2.0; print(gas, mH/gas); assert sig(mH/gas,2)==6.1
wt=5*1.0/(40.1+5*58.7+5*1.0)*100; print(wt); assert sig(wt,2)==1.5
print('OK')
