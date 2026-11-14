# INO-PPT-058 verification: two-wavelength simultaneous determination of Cu and Ni ammines
from sympy import symbols, solve, Rational as Q
from math import log10,floor
def sig(v,n): return round(v,-int(floor(log10(abs(v))))+n-1)
x,y=symbols('x y')
s=solve([52*x/10+Q(4)*y/10-Q('0.403'), 3*x+Q('6.5')*y-Q('0.386')],[x,y]); X,Y=float(s[x]),float(s[y])
print(X,Y); assert sig(X,3)==0.0756 and sig(Y,3)==0.0245
cu=X*0.1*63.5/0.8*100; ni=Y*0.1*58.7/0.8*100; zn=100-cu-ni; print(cu,ni,zn)
assert sig(cu,2)==60 and sig(ni,2)==18 and sig(zn,2)==22
A600=52*X+4*Y; print('A600 undiluted',A600); assert A600>1.0 and sig(A600,2)==4.0
r1=(4.03/52)/X; r2=(0.386/6.5)/Y; print(r1,r2); assert sig(r1,2)==1.0 and sig(r2,2)==2.4
f600=4*Y/A600; f360=3*X/0.386; print(f600,f360); assert round(f600*100)==2 and round(f360*100)==59
print('OK')
