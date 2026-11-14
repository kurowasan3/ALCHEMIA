# THE-EQL-159 verification: uncompetitive inhibition
from sympy import symbols, solve, simplify, Rational as Q
E,ES,ESI,S,I,Km,Ki,E0,k2=symbols('E ES ESI S I Km Ki E0 k2',positive=True)
sol=solve([E-Km*ES/S, ESI-ES*I/Ki, E0-E-ES-ESI],[E,ES,ESI],dict=True)[0]
v=simplify(k2*sol[ES]); a=1+I/Ki
assert simplify(v-(k2*E0*S/(Km+a*S)))==0
print('v =',v)
Vm,K,al=10,0.40,3
f=lambda S,alpha: Vm*S/(K+alpha*S)
data=[(0.20,f(0.20,1),f(0.20,al)),(4.0,f(4.0,1),f(4.0,al))]; print(data)
# from rounded data
d=[(0.20,3.33,2.00),(4.0,9.09,3.23)]
def line(i):
    x1,y1=1/d[0][0],1/d[0][i]; x2,y2=1/d[1][0],1/d[1][i]
    s=(y1-y2)/(x1-x2); return s, y2-s*x2
s0,b0=line(1); s1,b1=line(2); print(s0,b0,s1,b1)
Vmax=1/b0; Kmf=s0*Vmax; alpha=b1/b0; Kif=0.30/(alpha-1); print(Vmax,Kmf,alpha,Kif)
assert round(Vmax)==10 and round(Kmf,2)==0.40 and round(Kif,2)==0.15
for x in (0.1,1,10):
    u=(1+x)/(1+3*x); c=(1+x)/(3+x); print(x,round(u,2),round(c,2))
print('OK')
