# POL-PRT-076 verification: ESI multiply charged ions
from sympy import symbols, solve, Rational
mH=1.007; mNa=22.990
M=10846.0
peaks={z:(M+z*mH)/z for z in range(8,14)}
for z,p in peaks.items(): print(z,round(p,1))
p1=round(peaks[10],1); p2=round(peaks[11],1); p3=round(peaks[12],1)
print('given',p1,p2,p3)
z=(p2-mH)/(p1-p2); print('z from p1,p2',z)
zi=round(z); Mcalc=zi*(p1-mH); Mcalc2=(zi+1)*(p2-mH); Mcalc3=(zi+2)*(p3-mH)
print('M',Mcalc,Mcalc2,Mcalc3,'mean',(Mcalc+Mcalc2+Mcalc3)/3)
assert zi==10 and abs(Mcalc-M)<1.0
# symbolic derivation
Ms,zs,a,b=symbols('M z a b',positive=True)
sol=solve([Ms+zs*mH-zs*a, Ms+(zs+1)*mH-(zs+1)*b],[Ms,zs],dict=True)[0]
print('z formula',sol[zs])
# isotope spacing
print('spacing z=10',1.003/10,'z=11',1.003/11)
# peptide
Mp=1570.8
d2=(Mp+2*mH)/2; d3=(Mp+3*mH)/3; na=(Mp+mH+mNa)/2
print('pep 2+',d2,'3+',d3,'Na adduct 2+',na,'diff',na-d2)
assert round(d2,1)==786.4 and round(d3,1)==524.6 and round(na,1)==797.4
# back-calc from given
print('M from 786.4',2*(786.4-mH),'from 524.6',3*(524.6-mH))
# albumin
import math
zmin=66400/(2000-mH); print('zmin',zmin, math.ceil(zmin))
assert math.ceil(zmin)==34
print('OK')
