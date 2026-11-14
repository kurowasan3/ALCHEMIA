# THE-ELC-083 検証：EDI スタックの電流収支
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
print(check_equation('H2O -> H+ + OH-'), check_equation('H+ + OH- -> H2O'))
F=9.65e4
c=2.0e-3; Qd=0.50/60; Qc=0.050/60   # L/s
flux=c*Qd
print('(1) flux',flux,'Imin',F*flux)
I=2.0; e=I/F
print('(2) eff',flux/e,' water split mol/s',e-flux)
I3=1.2; print('(3) cout',c-I3/F/Qd)
print('(4) concentrate out',c+flux/Qc)
print('(5) total water split per 10 cells per hour mol', 10*(e-flux)*3600)
