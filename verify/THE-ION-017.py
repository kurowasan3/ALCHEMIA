import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from sympy import Rational as R, symbols, solve, log, N
a=symbols('a', positive=True); K2=R(1,100); C=R(3,1000)
al=solve(C*(1+a)*a-K2*(1-a), a); print(al); assert al==[R(2,3)]
h=C*(1+R(2,3)); print(h, N(-log(h,10),4), N(-log(2*C,10),4)); assert h==R(5,1000)
Ceq=solve(C_:=symbols('C')*R(3,2)-K2)  # 3/2 C = K2
print(Ceq, N(Ceq[0])); assert abs(float(Ceq[0])-6.67e-3)<1e-5
c=R(1,100); al2=solve(c*(1+a)*a-K2*(1-a), a)[0]; print('alpha at 0.01', N(al2)); assert al2<R(1,2)
assert not check_item('THE-ION-017')
