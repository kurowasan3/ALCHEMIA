import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from sympy import Rational as R
nAg=R(100,1000)*R(30,1000); nSCN=R(5,100)*R(20,1000); nX=nAg-nSCN
assert nX==R(2,1000); c=nX/R(20,1000); assert c==R(1,10)
M=R(119,100)/(nX*5); print('M=',M); assert M==119 and M-39==80
for X,m in (('Cl',35.5),('Br',80),('I',127)): print('K'+X, 39+m)
# 誤差の向き：Cl- が Ag+ を余分に消費 -> SCN 減 -> nX 計算値増
assert (nAg-(nSCN-R(1,10000))) > nX
assert not check_item('INO-PPT-019')
