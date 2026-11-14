# THE-BND-107 verification: spring model of diatomic IR absorption
from sympy import Rational as R, sqrt, nsimplify
c=3.00e10; pi=3.14; NA=6.02e23
m={'H':1.008,'D':2.014,'F':19.00,'35':34.97,'37':36.97,'Br':78.92}
mu=lambda a,b:a*b/(a+b)
def k(nu,a,b): return mu(a,b)/1000/NA*(2*pi*c*nu)**2
print('mu H35Cl kg', mu(m['H'],m['35'])/1000/NA)
kH=k(2886,m['H'],m['35']); kF=k(3962,m['H'],m['F']); kB=k(2559,m['H'],m['Br'])
print('k HCl,HF,HBr', kH,kF,kB)
assert round(kH)==481 and round(kF)==886 and abs(kB-384)<1
assert kF>kH>kB and 568>432>366
r=lambda a1,b1,a2,b2:(mu(a1,b1)/mu(a2,b2))**0.5
n37=2886*r(m['H'],m['35'],m['H'],m['37']); print('H37Cl',n37,'gap',2886-n37)
nD=2886*r(m['H'],m['35'],m['D'],m['35']); nDF=3962*r(m['H'],m['F'],m['D'],m['F'])
print('D35Cl',nD,'DF',nDF,'ratios',nD/2886,nDF/3962)
assert abs(n37-2883.8)<0.05 and round(nD)==2070 and round(nDF)==2873
print('mu range', mu(m['H'],m['F']), mu(m['H'],m['Br']))
print('OK')
