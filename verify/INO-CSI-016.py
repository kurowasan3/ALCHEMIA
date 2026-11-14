# INO-CSI-016 検証：Na 形 → NH4 形 → H 形ゼオライトの量的関係
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['NaAlSi2O6 + NH4+ -> NH4AlSi2O6 + Na+','NH4AlSi2O6 -> HAlSi2O6 + NH3','NaAlSi2O6 + 4HCl -> NaCl + AlCl3 + 2SiO2 + 2H2O']:
    assert check_equation(e) is None, e
H,N,O,Na,Al,Si=1.0,14,16,23,27,28
MNa=Na+Al+2*Si+6*O; MNH4=N+4*H+Al+2*Si+6*O; MH=H+Al+2*Si+6*O
n=10.0/MNa
print(MNa,MNH4,MH,n, n*22.4, n*MNH4, n*MH)
assert (MNa,MNH4,MH)==(202,197,180)
assert round(n*22.4,2)==1.11 and round(n*MNH4,2)==9.75 and round(n*MH,2)==8.91
# 酸点の密度の比較（1 g あたり）：Si/Al=2 と Si/Al=5（NaAlSi5O12, 脱水物）
MH5=H+Al+5*Si+12*O
print('mmol/g', 1000/MH, 1000/MH5)
print('OK')
