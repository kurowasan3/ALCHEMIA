import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
F = 9.65e4
ne = 0.448/22.4; I = ne*F/3860
nCl = 0.200*0.200; assert nCl > ne
nH = 0.500*2*30.0/1000; nK = nH - ne
VH2 = nK/2*22.4; VCa = (ne + 2*nK)/2/0.500*1000
print(ne, I, nK, nK*39.1, VH2, VCa)
assert round(I, 3) == 0.500 and round(nK*39.1, 3) == 0.391 and round(VH2, 3) == 0.112 and round(VCa, 1) == 40.0
for e in ['2KCl + 2H2O -> 2KOH + H2 + Cl2', '2H2O + 2e- -> H2 + 2OH-', '2Cl- -> Cl2 + 2e-', '2K + 2H2O -> 2KOH + H2', 'Ca + 2H2O -> Ca(OH)2 + H2']:
    assert check_equation(e) is None, e
print('OK')
