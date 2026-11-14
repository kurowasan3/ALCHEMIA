# INO-CSI-015 検証：前駆体高分子から得られる SiC の質量（反応式の釣り合いと数値）
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R
H,C,O,Si=R('1.0'),R('12.0'),R('16.0'),R('28.1')
eqs=['Si(CH3)2 -> SiC + CH4 + H2','SiH(CH3) -> SiC + 2H2','2Si(CH3)2 + 9O2 -> 2SiO2 + 4CO2 + 6H2O']
for e in eqs:
    assert check_equation(e) is None, e
P=Si+2*C+6*H; Q=Si+C+4*H; SiC=Si+C; SiO2=Si+2*O
print('式量 P,Q,SiC,SiO2', P,Q,SiC,SiO2)
yP=100*SiC/P; yQ=100*SiC/Q; ox=100*SiO2/P
print(float(yP), float(yQ), float(ox))
assert round(float(yP),1)==69.0 and round(float(yQ),1)==90.9 and round(float(ox))==103
print('OK')
