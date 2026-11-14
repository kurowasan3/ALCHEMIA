import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
from fractions import Fraction as F
Qf = lambda n, h, a: F(a)**2/(F(n)*F(h)**3)
K = Qf('0.5', '0.5', '1.0'); print('K', K); assert K == 16
assert Qf('1.0', '1.0', '2.0') == 4 and Qf('1.0', '1.0', '2.0') < K
assert Qf('0.5', '0.5', '1.0') == K        # 2.0 L に 1.0,1.0,2.0 mol
assert Qf('0.5', '0.5', '2.0') == 64 and 64 > K
assert check_item('THE-EQL-024')
