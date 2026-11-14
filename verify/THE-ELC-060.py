# THE-ELC-060 検証：LiCn と水の反応による組成決定・充電の割合
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
Li, C = 6.94, 12.0
nLi = 2*0.284/22.4; M = 2.00/nLi; n = (M-Li)/C
print('nLi', nLi, 'M', M, 'n', n)
assert round(n) == 6 and abs(n-6) < 0.05
print('HCl mL', nLi/0.500*1000)
assert round(nLi/0.500*1000, 1) == 50.7
nLi2 = 2*0.142/22.4; mLiC6 = nLi2*(Li+6*C); mC_free = 2.00 - mLiC6
Ctot = nLi2*6 + mC_free/C; soc = nLi2/(Ctot/6)
print('LiC6 g', mLiC6, 'free C g', mC_free, 'SOC', soc)
assert round(soc*100) == 48
for e in ['2LiC6 + 2H2O -> 2LiOH + 12C + H2', 'LiOH + HCl -> LiCl + H2O', '6C + Li+ + e- -> LiC6']:
    assert check_equation(e) is None, e
print('OK')
