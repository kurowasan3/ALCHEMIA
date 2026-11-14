# POL-PRT-055 検証：側鎖アミドの加水分解とアンモニアの定量
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q
from rdkit import Chem
assert check_equation('CH3CONH2 + H2O + H+ -> CH3COOH + NH4+') is None
Hp = Q('0.100') * 2 * Q('50.0'); OH = Q('0.100') * Q('60.0')
nh3 = Hp - OH; print(nh3, nh3 / 2)
assert nh3 == 4 and nh3 / 2 == 2
# Ala-Lys-Asn-Gln（配列は一例）
pep = Chem.MolFromSmiles('C[C@H](N)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCC(N)=O)C(=O)O')
nN = sum(a.GetSymbol() == 'N' for a in pep.GetAtoms()); print('N atoms', nN); assert nN == 7
charge = +1 - 1 + 1  # N末端, C末端, Lys側鎖（pH 6.0）
assert charge == 1
pI = {'Ala': 6.0, 'Lys': 9.7, 'Asp': 2.8, 'Glu': 3.2}
anode = sorted(k for k, v in pI.items() if v < 6.0 - 0.5); print(anode); assert anode == ['Asp', 'Glu']
print('OK')
