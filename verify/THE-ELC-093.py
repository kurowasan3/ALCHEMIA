# THE-ELC-093 検証：現像の局部電池モデル（AgBr とヒドロキノン）
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
eqs = ['AgBr + e- -> Ag + Br-',
       'C6H6O2 -> C6H4O2 + 2H+ + 2e-',
       '2AgBr + C6H6O2 -> 2Ag + C6H4O2 + 2H+ + 2Br-',
       'AgBr + 2S2O3^2- -> [Ag(S2O3)2]^3- + Br-',
       '2AgBr -> 2Ag + Br2']
for e in eqs:
    r = check_equation(e); print(e, r); assert r is None
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
print(CalcMolFormula(Chem.MolFromSmiles('Oc1ccc(O)cc1')), CalcMolFormula(Chem.MolFromSmiles('O=C1C=CC(=O)C=C1')))
nAg = 0.216/107.9
nHQ = nAg/2
mHQ = nHQ*(6*12.0 + 6*1.0 + 2*16.0)
Q = nAg*9.65e4
print('nAg', nAg, 'HQ mol', nHQ, 'HQ g', mHQ, 'Q', Q)
assert f'{mHQ:.3f}' == '0.110' and f'{Q:.3g}' == '193'
