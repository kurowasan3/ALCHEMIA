import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D
H, C, O = 1.0, 12, 16
def M(f): return D.CalcMolFormula(Chem.MolFromSmiles(f))
assert M('C=CC(=O)OC(C)=C') == 'C6H8O2' and M('C=CC(=O)OCC=C') == 'C6H8O2'
Mest = 6*C + 8*H + 2*O; Macid = 3*C + 4*H + 2*O
n0 = 11.2/Mest; nH2 = 0.340/(2*H)
a = nH2 - n0          # 2a + (n0 - a) = nH2
hyd = n0 - a
print(n0, nH2, a, hyd, hyd/n0, hyd*Macid)
assert abs(hyd/n0 - 0.30) < 1e-9 and round(hyd*Macid, 2) == 2.16
assert check_equation('CH2=CHCOOC(CH3)=CH2 + H2O -> CH2=CHCOOH + CH3COCH3') is None
print('OK')
