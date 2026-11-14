# ORG-CBX-046 検証：水素の吸収量と還元的オゾン分解による不飽和脂肪酸の構造
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import MolWt
X = 'CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)O'          # 9,12-オクタデカジエン酸
prods = {'E': 'CCCCCCO', 'F': 'OCCCO', 'G': 'OCCCCCCCCC(=O)O'}
fx = CalcMolFormula(Chem.MolFromSmiles(X)); print('X', fx); assert fx == 'C18H32O2'
M = 18*12.0 + 32*1.0 + 2*16.0; print('M(X) =', M); assert M == 280.0
nX = 0.280/M; nH2 = 44.8e-3/22.4; print('H2/X =', nH2/nX); assert round(nH2/nX, 6) == 2
for k, s in prods.items(): print(k, CalcMolFormula(Chem.MolFromSmiles(s)))
assert [CalcMolFormula(Chem.MolFromSmiles(s)) for s in prods.values()] == ['C6H14O', 'C3H8O2', 'C9H18O3']
# 炭素数の保存 6 + 3 + 9 = 18
assert 6 + 3 + 9 == 18
# 還元的オゾン分解を機械的に行う（C=C を切って両側を CH2OH にする）
def reductive_ozonolysis(smi):
    m = Chem.RWMol(Chem.MolFromSmiles(smi)); Chem.Kekulize(m, clearAromaticFlags=True)
    dbl = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()
           if b.GetBondType() == Chem.BondType.DOUBLE and m.GetAtomWithIdx(b.GetBeginAtomIdx()).GetSymbol() == 'C'
           and m.GetAtomWithIdx(b.GetEndAtomIdx()).GetSymbol() == 'C']
    for a, b in dbl:
        m.RemoveBond(a, b)
        for at in (a, b):
            o = m.AddAtom(Chem.Atom(8)); m.AddBond(at, o, Chem.BondType.SINGLE)
    mol = m.GetMol(); Chem.SanitizeMol(mol)
    frags = Chem.GetMolFrags(mol, asMols=True)
    return sorted(Chem.MolToSmiles(f) for f in frags)
r1 = reductive_ozonolysis(X); print('X ->', r1)
assert sorted(r1) == sorted(Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in prods.values())
L = 'CC/C=C\\C/C=C\\C/C=C\\CCCCCCCC(=O)O'      # 9,12,15-オクタデカトリエン酸（リノレン酸）
r2 = reductive_ozonolysis(L); print('linolenic ->', r2)
assert r2.count('OCCCO') == 2 and 'CCCO' in r2 and len(set(r2)) == 3
# 三重結合1個（例：9-オクタデシン酸）なら切断箇所1つで断片は2つ（ジオールは生じない）
print('OK')
