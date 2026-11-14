import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from math import log10
pP=-log10(1e-9); pQ=-log10(2e-5); print(pP,pQ)
rng={'P':(pP-1,pP+1),'Q':(round(pQ,1)-1,round(pQ,1)+1)}; print(rng)
assert rng['P']==(8.0,10.0) and abs(rng['Q'][0]-3.7)<1e-9
assert abs(1e-7/1e-9-100)<1e-9
inside=lambda r,p: r[0]<=p<=r[1]
assert inside(rng['P'],8.7) and not inside(rng['Q'],8.7)
assert inside(rng['Q'],5.3) and not inside(rng['P'],5.3)
assert all(3<=r[0] and r[1]<=11 for r in rng.values())
assert not check_item('THE-ION-014')
