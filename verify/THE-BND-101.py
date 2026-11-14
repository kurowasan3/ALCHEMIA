# THE-BND-101 verification: electrons around the central atom, formal charges, equations
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
val = {'P':5,'S':6,'Xe':8,'Cl':7,'Si':4,'N':5,'O':6,'C':4,'I':7}
def around(center, nbonds, charge=0):
    e = val[center] - charge          # anion charge adds electrons
    lone = e - nbonds                 # all bonds single, one electron from center each
    assert lone >= 0 and lone % 2 == 0
    return 2*nbonds + lone
res = {'PCl5':around('P',5),'SF4':around('S',4),'XeF2':around('Xe',2),'ClF3':around('Cl',3),
       '[SiF6]2-':around('Si',6,-2),'NF3':around('N',3)}
print(res)
assert res == {'PCl5':10,'SF4':10,'XeF2':10,'ClF3':10,'[SiF6]2-':12,'NF3':8}
# (2) second-period centers exceeding 8
for f,c,n in [('PF5','P',5),('NF5','N',5),('SF6','S',6),('OF4','O',4),('IF7','I',7)]:
    print(f, around(c,n))
assert around('N',5) == 10 and around('O',4) == 10
# (3) H2SO4 all single bonds: S 4 bonds, 0 lone; terminal O 1 bond + 3 lone pairs
fc = lambda v, lone_e, bonds: v - lone_e - bonds
assert fc(6,0,4) == 2 and fc(6,6,1) == -1 and fc(6,4,2) == 0
assert 2*4 == 8
# (5) CF4 already 8
assert around('C',4) == 8
for e in ['SiF4 + 2F- -> [SiF6]^2-', 'SiO2 + 6HF -> H2SiF6 + 2H2O']:
    r = check_equation(e); print(e, r); assert r is None
print('OK')
