import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from math import sqrt
assert abs(10.0+26.5+1.0+12.5-50.0)<1e-9
c=0.0200*12.5/10.0; print(c); assert abs(c-0.025)<1e-12
cro4=0.10*1.0/50.0; assert abs(cro4-2.0e-3)<1e-12
a1=1.8e-10/5.0e-3; a2=sqrt(2.0e-12/2.0e-3); print(a1,a2); assert round(a1,10)==3.6e-8 and abs(a2-3.16e-5)<1e-7
cl=1.8e-10/3.16e-5; print(cl); assert f'{cl:.1e}'=='5.7e-06'
print('excess Ag mL', a2*0.050/0.0200*1000)
assert not check_item('THE-ION-016')
