# ORG-ARO-072 検証：燃焼データ・水素の吸収量・候補構造
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import MolWt
H, C, O = 1.0, 12.0, 16.0
mC = 19.8*C/(C+2*O); mH = 8.10*2*H/(2*H+O)
print('C', mC, 'H', mH, 'sum', mC+mH, 'molar C:H', (mC/C)/(mH/H))
assert abs(mC+mH-6.30) < 1e-9 and abs((mC/C)/(mH/H)-0.5) < 1e-9
n = max(k for k in range(1, 20) if 14*k <= 130); print('n max', n)
nA = 0.448/22.4/2; MA = 1.16/nA; print('A molar mass', MA)
assert abs(MA - (9*C+8*H)) < 1e-9
cands = {'3-phenylpropyne': 'C#CCc1ccccc1', '1-phenylpropyne': 'CC#Cc1ccccc1',
         'o-ethynyltoluene': 'C#Cc1ccccc1C', 'm-': 'C#Cc1cccc(C)c1', 'p-': 'C#Cc1ccc(C)cc1'}
for k, s in cands.items():
    m = Chem.MolFromSmiles(s); print(k, CalcMolFormula(m))
    assert CalcMolFormula(m) == 'C9H8'
# 末端アルキン（≡C–H）を持つか
patt = Chem.MolFromSmarts('[CH1]#C')
term = [k for k, s in cands.items() if Chem.MolFromSmiles(s).HasSubstructMatch(patt)]
print('terminal alkyne:', term)
# 一置換ベンゼン（KMnO4 で安息香酸）
mono = [k for k in term if k == '3-phenylpropyne']
assert term == ['3-phenylpropyne', 'o-ethynyltoluene', 'm-', 'p-']
for e in ['2C9H18 + 27O2 -> 18CO2 + 18H2O', 'C6H5CH2C#CH + 2H2 -> C6H5CH2CH2CH3',
          'C6H5CH2CH2CH3 + 3H2 -> C6H11CH2CH2CH3']:
    assert check_equation(e) is None, (e, check_equation(e))
print('D=C9H18 (propylcyclohexane):', CalcMolFormula(Chem.MolFromSmiles('CCCC1CCCCC1')))
print('OK')
