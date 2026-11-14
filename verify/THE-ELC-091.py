# THE-ELC-091 verification: graphite vs silicon anode capacity
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
F=9.65e4; Li=6.94; C=12.0; Si=28.1
g1=F/(6*C)/3.6; s1=4.4*F/Si/3.6; print(g1,s1,s1/g1)
assert round(g1)==372 and round(s1,-1)==4200 and round(s1/g1)==11
g2=F/(6*C+Li)/3.6; s2=4.4*F/(Si+4.4*Li)/3.6; print(g2,s2,s2/g2)
assert round(g2)==340 and round(s2,-1)==2010 and abs(s2/g2-5.9)<0.05
Vg=6*C/2.25*1.10; Vs=Si/2.33*4.0
g3=F/3.6/Vg; s3=4.4*F/3.6/Vs; print(Vg,Vs,g3,s3,s3/g3)
assert round(g3)==762 and round(s3,-1)==2440 and abs(s3/g3-3.2)<0.05
for e in ['6C + Li+ + e- -> LiC6','5Si + 22Li+ + 22e- -> Li22Si5']:
    assert check_equation(e) is None
print('OK')
