# THE-THC-072 verification: photodissociation threshold and chain reaction quantum yield
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
h=6.63e-34; c=3.00e8; NA=6.02e23
lCl=h*c*NA/243e3*1e9; lH=h*c*NA/436e3*1e9
print('lambda Cl2,H2 nm', lCl, lH); assert round(lCl)==493 and round(lH)==275
assert 600>lCl and 450<lCl
Eph=h*c/450e-9; nph=2.0e-3*100/Eph/NA
print('Eph',Eph,'photon mol',nph,'QY',1.5e-2/nph)
assert abs(nph-7.5e-7)<0.05e-7 and round(1.5e-2/nph,-3)==20000
print('steps', 436-432, 243-432, 2*432-436-243)
q=1.5e-2/2*(2*432-436-243); print('heat kJ',q,'ratio',q*1e3/0.20); assert round(q,1)==1.4
for e in ['Cl2 -> 2Cl','Cl + H2 -> HCl + H','H + Cl2 -> HCl + Cl','H2 + Cl2 -> 2HCl']:
    assert check_equation(e) is None, e
print('OK')
