import sys; sys.path.insert(0, '/home/claude/chem/work/g04')
from vlib import check_item, approx
approx(0.0500*0.200*126, 1.26); approx(2*0.0500*10.0/8.00, 0.125)
c = lambda acid_mol_factor, vread: 2*0.0500*10.0*acid_mol_factor/vread
true = c(1, 8.00)
assert c(1, 8.00*0.95) > true      # (a) 酸が薄まり滴下量減
assert c(1, 8.00) == true          # (b)
assert c(1, 8.00*1.05) < true      # (c) NaOH が薄まり滴下量増
assert c(1, 8.00+0.05) < true      # (d) 気泡分だけ読みが大
assert check_item('THE-ACB-013')
