# THE-EQL-148 verification: irreversible vs competitive inhibition
from sympy import Rational as Q
Vm=Q(1,2); Km=2
v=lambda V,K,S: V*S/(K+S)
VX=Vm*Q(6,10); KY=Km*(1+Q(3,10)/Q(1,10))
for S in (2,50):
    a,b,c=v(Vm,Km,S),v(VX,Km,S),v(Vm,KY,S)
    print(S,float(a),float(b),float(c),float(b/a),float(c/a))
assert VX==Q(3,10) and KY==8
assert abs(float(v(VX,Km,2))-0.15)<1e-12 and abs(float(v(Vm,KY,2))-0.10)<1e-12
assert round(float(v(Vm,KY,50)),2)==0.43 and round(float(v(VX,Km,50)),2)==0.29 and round(float(v(Vm,Km,50)),2)==0.48
print('OK')
