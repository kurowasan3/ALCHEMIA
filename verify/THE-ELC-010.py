import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from math import log10
E={'A':0.0}; E['B']=E['A']+0.78; E['C']=E['A']-0.50; E['D']=E['A']+1.24
order=sorted(E, key=E.get); print(order); assert order==['C','A','B','D']
assert abs(E['D']-E['C']-1.74)<1e-9 and abs(E['D']-E['B']-0.46)<1e-9
dil=0.059/2*log10(0.010); conc=0.0; assert dil<conc; print('conc cell emf', conc-dil)
assert not check_item('THE-ELC-010')
