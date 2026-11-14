# INO-TRN-047 verification: chlorine leaching of Ni3S2/Cu2S matte, Cl2 recycling, Mond process
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, solve
for e in ['Ni3S2 + 3Cl2 -> 3NiCl2 + 2S','Cu2S + 2Cl2 -> 2CuCl2 + S','NiCl2 -> Ni + Cl2','Ni + 4CO -> Ni(CO)4']:
    assert check_equation(e) is None, e
x,y=symbols('x y')
s=solve([240.1*x+159.0*y-2000, 32.0*(2*x+y)-480],[x,y]); X,Y=float(s[x]),float(s[y])
def sig(v,n): 
    from math import log10,floor; return round(v,-int(floor(log10(abs(v))))+n-1)
print(X,Y)
assert sig(X*240.1/1000,2)==1.2 and sig(Y*159.0/1000,2)==0.81 and sig(3*X*58.7/1000,2)==0.87
assert sig(3*X+2*Y,2)==25 and sig(2*Y,2)==10
assert sig(4*3*X*22.4,2)==1300
# Mond: Delta H from formation enthalpies (kJ/mol): Ni(CO)4(g) -602.9, CO(g) -110.5
print('dH',-602.9-4*(-110.5)); assert abs((-602.9-4*(-110.5))-(-161))<1
print('OK')
