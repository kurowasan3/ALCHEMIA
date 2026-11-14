import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
from sympy import Rational as R, symbols, solve, sqrt, N
x=symbols('x', positive=True)
print(N(sqrt(R(2,10**10)))); assert abs(float(sqrt(2e-10))-1.41e-5)<1e-7
Kf=R(125,10)*10**6; K=R(2,10**10)*Kf; assert K==R(25,10**4)
xs=[s for s in solve(x**2-K*(1-2*x)**2, x) if s<R(1,2)][0]; print(N(xs))
K2=R(5,10**13)*Kf; ys=[s for s in solve(x**2-K2*(1-2*x)**2, x) if s<R(1,2)][0]; print(N(ys), N(xs/ys))
assert f'{float(xs):.1e}'=='4.5e-02' and f'{float(ys):.1e}'=='2.5e-03'
assert not check_item('THE-ION-019')
