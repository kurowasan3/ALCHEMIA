# INO-PPT-054 ジメチルグリオキシムによるニッケルの重量分析
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
H,C,N,O,Ni,Cu = 1.0,12.0,14.0,16.0,58.7,63.5
dmgH2 = Chem.MolFromSmiles('CC(=NO)C(C)=NO'); print('dmgH2', CalcMolFormula(dmgH2))   # C4H8N2O2
M_complex = 8*C + 14*H + 4*N + 4*O + Ni    # Ni(C4H7N2O2)2
print('formula mass Ni(dmgH)2 =', M_complex)
for eq in ['Ni^2+ + 2C4H8N2O2 -> Ni(C4H7N2O2)2 + 2H+', 'H+ + NH3 -> NH4+',
           'Ni^2+ + 2C4H8N2O2 + 2NH3 -> Ni(C4H7N2O2)2 + 2NH4+']:
    print(eq, check_equation(eq))
sample = 0.250; ppt = 0.302
mNi = ppt*Ni/M_complex; print('Ni mass', mNi, 'Ni % =', mNi/sample*100, 'Cu % =', 100-mNi/sample*100)
print('ppt per mg Ni =', M_complex/Ni)
# 必要な dmgH2 の最小量
print('dmgH2 min (g) =', 2*mNi/Ni*(4*C+8*H+2*N+2*O))
