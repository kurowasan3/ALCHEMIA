# THE-EQL-161 verification: enantiomer ratio of a chiral pollutant in a well-mixed lake (steady state and decay)
from sympy import symbols, Function, dsolve, Eq, exp, simplify, solve, log, N
t=symbols('t',positive=True); F,V,ko,kd=symbols('F V k_o k_d',positive=True)
P=Function('P'); Mn=Function('M')
# d[+]/dt = F/(2V) - (ko+kd)[+] ; d[-]/dt = F/(2V) - ko[-]
Pss=solve(Eq(F/(2*V)-(ko+kd)*symbols('x'),0),symbols('x'))[0]
Mss=solve(Eq(F/(2*V)-ko*symbols('y'),0),symbols('y'))[0]
r=simplify(Pss/Mss); print('steady ratio',r); assert simplify(r-ko/(ko+kd))==0
# after stop
solP=dsolve(Eq(P(t).diff(t),-(ko+kd)*P(t)),P(t),ics={P(0):Pss}).rhs
solM=dsolve(Eq(Mn(t).diff(t),-ko*Mn(t)),Mn(t),ics={Mn(0):Mss}).rhs
rt=simplify(solP/solM); print('ratio(t)',rt); assert simplify(rt-r*exp(-kd*t))==0
# numbers lake A
koA=0.20; rA=0.80; kdA=koA*(1/rA-1); print('kd',kdA); assert abs(kdA-0.050)<1e-12
# lake B continuous input ko=0.020
rB=0.020/(0.020+kdA); print('rB',rB); assert round(rB,2)==0.29
# time to reach 0.50 after stop
import math
tA=math.log(0.80/0.50)/kdA; print('t',tA,'ln1.6',math.log(1.6)); assert round(tA,1)==9.4
# remaining total fraction
Mrem=math.exp(-koA*tA); Prem=0.8*math.exp(-(koA+kdA)*tA)
tot=(Mrem+Prem)/1.8; print('Mrem',Mrem,'Prem',Prem,'total fraction',tot); assert round(tot,2)==0.13
# lake B after stop: same ratio decay rate kd
tB=math.log(rB/(rB*0.625))/kdA; print('tB to fall to 62.5% of its value',tB)
# with ln values supplied: ln1.6=0.47
print('t using ln1.6=0.47', 0.47/0.05)
print('OK')
