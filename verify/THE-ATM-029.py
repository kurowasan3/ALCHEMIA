# THE-ATM-029 verification: mass defect of D-T fusion vs combustion of hydrogen
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
u=1.661e-27; c=3.00e8; NA=6.02e23
D,T,He,n=2.01410,3.01605,4.00260,1.00866
dm=D+T-He-n
print('dm u',dm,'fraction %',dm/(D+T)*100)
assert abs(dm-0.01889)<1e-6
dm_kg=dm*u; E1=dm_kg*c**2
print('E per reaction J',E1)
Emol=dm*1e-3*c**2   # kg/mol * c^2
print('E per mol J',Emol,'per g fuel J',Emol/(D+T))
assert round(Emol/1e12,2)==1.70 and round(Emol/(D+T)/1e11,2)==3.38
# combustion 2H2+O2->2H2O, 572 kJ ; per g of mixture (4.0+32)=36 g
Ecomb_g=572e3/36.0
print('comb per g J',Ecomb_g,'ratio',Emol/(D+T)/Ecomb_g)
assert round(Ecomb_g/1e4,2)==1.59 and round(Emol/(D+T)/Ecomb_g/1e7,1)==2.1
dm_comb=572e3/c**2*1e3  # g
print('mass loss on combustion of 2 mol H2 (g)',dm_comb)
assert round(dm_comb*1e9,1)==6.4
# same per 1 g of mixture
print('fraction of mass lost combustion',572e3/c**2/0.036)
assert check_equation('2H2 + O2 -> 2H2O') is None
print('OK')
