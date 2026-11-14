import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
import math
Kw = 1e-14
K1c, K2c = 4.5e-7, 4.7e-11; K1p, K2p, K3p = 7.1e-3, 6.3e-8, 4.5e-13
approx(Kw/K1c, 2.2e-8, 0.02); assert Kw/K1c > K2c
approx(Kw/K1p, 1.4e-12, 0.02); assert Kw/K1p < K2p
approx(Kw/K2p, 1.6e-7, 0.02); assert Kw/K2p > K3p
approx(Kw/1.8e-5, 5.6e-10, 0.01)
pk = lambda k: -math.log10(k)
print('pH NaHCO3', (pk(K1c)+pk(K2c))/2, 'NaH2PO4', (pk(K1p)+pk(K2p))/2, 'Na2HPO4', (pk(K2p)+pk(K3p))/2)
assert check_item('THE-ION-010')
