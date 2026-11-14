# THE-BND-108 verification: iron allotropic transformations
from math import sqrt
NA=6.02e23; Fe=55.85
af=0.2904*sqrt(3)/sqrt(2); print('model a_f',af); assert round(af,3)==0.356
vb=(4/sqrt(3))**3/2; vf=(2*sqrt(2))**3/4; print('model ratio',vf/vb, 1-vf/vb); assert round((1-vf/vb)*100,1)==8.1
V={}
for n,a,z in [('b912',0.2904,2),('f912',0.3647,4),('f1394',0.3687,4),('b1394',0.2932,2)]:
    V[n]=a**3/z; rho=Fe/(NA*V[n]*1e-21); r=sqrt(3)*a/4 if z==2 else sqrt(2)*a/4
    print(n,V[n],rho,r)
d912=1-V['f912']/V['b912']; d1394=V['b1394']/V['f1394']-1
print(d912,d1394); assert round(d912*100,2)==0.97 and round(d1394*100,2)==0.58
rb=sqrt(3)*0.2904/4; rf=sqrt(2)*0.3647/4; print(rf/rb,(rf/rb)**3)
assert round(rb,3)==0.126 and round(rf,3)==0.129
print('OK')
