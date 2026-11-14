import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
import math
c1 = 0.020*20.0/10.0; c2 = 0.020*10.0/10.0; c3 = 0.020*25.0/10.0
approx(c1, 0.040); approx(c2, 0.020); approx(c3, 0.050)
approx(1e-3/c1, 0.025); approx(1e-3/c3, 0.020)
approx(14+math.log10(c2), 12.3, 0.003)           # P 初め pH
Ka = 2.5e-5; x = (-Ka+math.sqrt(Ka*Ka+4*Ka*c1))/2; print('R 初め pH', -math.log10(x)); approx(-math.log10(x), 3.0, 0.01)
cs = c1*10/30; oh = math.sqrt(1e-14/Ka*cs); print('R 中和点 pH', 14+math.log10(oh)); approx(14+math.log10(oh), 8.4, 0.01)
Kb = c3*0.02**2/(1-0.02); cn = c3*10/35; h = math.sqrt(1e-14/Kb*cn); print('Q 中和点 pH', -math.log10(h)); approx(-math.log10(h), 5.6, 0.01)
approx(-math.log10(0.020), 1.7, 0.01)
assert check_item('THE-ACB-011')
