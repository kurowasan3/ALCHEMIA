import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from itertools import permutations
rank={'Mg':0,'Al':1,'Fe':2,'Pb':3,'Cu':4,'Ag':5}
dil_HCl={'Mg','Al','Fe'}; hot_water={'Mg'}; passive={'Al','Fe'}; amph={'Al'}; dil_HNO3={'Mg','Al','Fe','Pb','Cu','Ag'}
sols=[]
for p in permutations(rank):
    P,Q,R,T,U,V=p
    ok = {P,Q,R}==dil_HCl and Q in hot_water and P in passive and R in passive and P in amph and R not in amph \
         and {T,U,V}<=dil_HNO3 and rank[U]>rank[V] and rank[V]>rank[T]
    if ok: sols.append(p)
print(sols); assert sols==[('Al','Mg','Fe','Pb','Ag','Cu')]
from chemlib import check_equation
assert check_equation('Pb + 2HCl -> PbCl2 + H2') is None
assert not check_item('THE-RDX-015')
