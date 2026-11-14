# ORG-HYD-047 検証：環状トリエンのオゾン分解と質量増加
from sympy import symbols, Rational
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import Descriptors
# C12Hh, k 個の C=C：8.10(144+h)^-1 * 32k = 4.80
sols = []
for k in range(1, 8):
    h = Rational(810*32*k, 480) - 144
    if h > 0 and h == int(h) and h % 2 == 0 and h <= 26:
        sols.append((k, int(h)))
print('solutions (k,h)', sols)
assert sols == [(3, 18)]
B = Chem.MolFromSmiles('C1=CCCC=CCCC=CCC1')
print(CalcMolFormula(B), [len(r) for r in B.GetRingInfo().AtomRings()])
assert CalcMolFormula(B) == 'C12H18'
# オゾン分解（還元的）：すべての C=C を切って両側を C=O にする
def ozonolysis(mol):
    rw = Chem.RWMol(mol)
    pairs = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
             if b.GetBondType() == Chem.BondType.DOUBLE and b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C']
    for a, b in pairs:
        rw.RemoveBond(a, b)
        for x in (a, b):
            o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(x, o, Chem.BondType.DOUBLE)
    m = rw.GetMol(); Chem.SanitizeMol(m)
    return m
m = ozonolysis(B)
prods = Chem.MolToSmiles(m).split('.')
print(prods)
assert set(prods) == {'O=CCCC=O'} and len(prods) == 3
print('succinaldehyde', CalcMolFormula(Chem.MolFromSmiles('O=CCCC=O')))
print('mass B', 8.10, 'n', 8.10/162, 'prod mass', 8.10/162*3*86, 'H2 mol', 8.10/162*3)
assert abs(8.10/162*3*86 - 12.90) < 1e-9
print('OK')
