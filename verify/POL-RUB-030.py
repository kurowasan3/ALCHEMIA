# POL-RUB-030 検証：ポリアセチレンの生成・ドーピングの反応式と構造
import sys; sys.path.insert(0,'/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['3I2 + 2e- -> 2I3-','Na -> Na+ + e-','CaC2 + 2H2O -> Ca(OH)2 + C2H2','3C2H2 -> C6H6','I2 + I- -> I3-','2C2H2 + 5O2 -> 4CO2 + 2H2O']:
    print(e, check_equation(e))
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
# オリゴマー（4単位）で共役（単結合・二重結合の交互）を確認
m=Chem.MolFromSmiles('C=CC=CC=CC=C'); print(CalcMolFormula(m))
# 26.0 g のポリアセチレン中の C=C の物質量
print('C=C mol', 26.0/26.0)
