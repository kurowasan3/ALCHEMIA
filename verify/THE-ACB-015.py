import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from math import log10
for c in (0.100, 0.0100):
    for V, exact in ((0.0400, None), (0.03995, 'before'), (0.04005, 'after')):
        pass
    before = -log10(c*0.050e-3/0.0400); after = 14 + log10(c*0.050e-3/0.0400)
    before_ex = -log10(c*0.050e-3/0.03995); after_ex = 14 + log10(c*0.050e-3/0.04005)
    print(c, round(before,2), round(after,2), round(before_ex,2), round(after_ex,2))
    assert round(before,1)==round(before_ex,1) and round(after,1)==round(after_ex,1)
assert round(-log10(1.25e-4),1)==3.9 and round(-log10(1.25e-5),1)==4.9
v = 20.0 - 4.0e-5*0.0400/0.0100*1000; print('V(MO)=', v); assert round(v,1)==19.8
print('pH60', -log10(1e-13**0.5)); assert abs(-log10(1e-13**0.5)-6.5)<1e-9
assert not check_item('THE-ACB-015')
