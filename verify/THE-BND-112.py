# THE-BND-112 verification: corundum/ilmenite stacking period, Z, electrostatic valence, density
from math import gcd, log10, floor
from fractions import Fraction as F
def lcm(*a):
    r=1
    for x in a: r=r*x//gcd(r,x)
    return r
def sig(v,n): return round(v,-int(floor(log10(abs(v))))+n-1)
assert F(6)*F(2,3)==4 and 4*F(3,6)==2
assert 2*F(2,6)+2*F(4,6)==2 and 4*F(2,6)==F(4,3) and 4*F(4,6)==F(8,3)
assert lcm(2,3)==6 and lcm(2,3,2)==6 and lcm(2,3,4)==12
O=6*3; M=F(6*3*2,3); assert (O,M)==(18,12) and M/2==6
NA=6.02e23
def rho(Mf,a,c): return 6*Mf/NA/(1.73/2*a*a*c)
r1=rho(2*55.8+3*16.0,5.04e-8,1.38e-7); r2=rho(55.8+47.9+3*16.0,5.09e-8,1.41e-7)
print(r1,r2); assert sig(r1,2)==5.2 and sig(r2,2)==4.8
# stacking sequence check: pairs (O layer, hole layer) repeat after 6
seq=[('AB'[i%2],'XYZ'[i%3]) for i in range(24)]
per=next(p for p in range(1,24) if all(seq[i]==seq[i+p] for i in range(24-p))); print(per, ''.join(a+b for a,b in seq[:6])); assert per==6
print('OK')
