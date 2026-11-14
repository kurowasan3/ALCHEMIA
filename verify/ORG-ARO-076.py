# ORG-ARO-076 検証：アセチル・ベンゾイルの 2 つのアシル基をもつイミド
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
D = {'CC(=O)N(C(=O)c1ccccc1)c1ccccc1': 'C15H13NO2', 'O=C(Nc1ccccc1)c1ccccc1': 'C13H11NO',
     'CC(=O)Nc1ccccc1': 'C8H9NO', 'CC(=O)N(C(C)=O)c1ccccc1': 'C10H11NO2'}
for s, f in D.items(): assert F(Chem.MolFromSmiles(s)) == f
for e in ['C15H13NO2 + H2O -> C13H11NO + C2H4O2', 'C15H13NO2 + H2O -> C8H9NO + C7H6O2',
          'C13H11NO + H2O -> C7H6O2 + C6H7N', 'C8H9NO + H2O -> C2H4O2 + C6H7N']:
    assert check_equation(e) is None
M = lambda c, h, n, o: 12*c + h + 14*n + 16*o
assert (M(15, 13, 1, 2), M(8, 9, 1, 1), M(13, 11, 1, 1)) == (239, 135, 197)
nX = Q('2.39')/239; nA = Q('0.54')/135; nB = nX - nA
print(nX, nA, float(nB*197), float(nA/nX*100))
assert round(float(nB*197), 1) == 1.2 and nA/nX == Q(2, 5)
# Y は 2 つの C-N 結合が等価：生成するアミドの種類
y = Chem.MolFromSmiles('CC(=O)N(C(C)=O)c1ccccc1'); r = list(Chem.CanonicalRankAtoms(y, breakTies=False))
carb = [a.GetIdx() for a in y.GetAtoms() if a.GetSymbol() == 'C' and any(n.GetSymbol() == 'N' for n in a.GetNeighbors()) and not a.GetIsAromatic()]
assert len({r[i] for i in carb}) == 1
print('OK')
