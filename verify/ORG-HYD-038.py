"""ORG-HYD-038 検証：マルコフニコフ水和 → ザイツェフ脱水"""
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

alk = {'a': 'C=C(C)CCC', 'b': 'CC(C)=CCC', 'c': 'C=CCC(C)C', 'd': 'C=C(C)C(C)C', 'e': 'C=CC(C)(C)C'}

def hydrate(smi):
    m = Chem.MolFromSmiles(smi)
    b = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE][0]
    a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
    oh = a2 if a1.GetTotalNumHs() > a2.GetTotalNumHs() else a1
    assert a1.GetTotalNumHs() != a2.GetTotalNumHs()
    rw = Chem.RWMol(m); rw.GetBondBetweenAtoms(a1.GetIdx(), a2.GetIdx()).SetBondType(Chem.BondType.SINGLE)
    o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(oh.GetIdx(), o, Chem.BondType.SINGLE)
    mm = rw.GetMol(); Chem.SanitizeMol(mm)
    return Chem.MolToSmiles(mm)

def dehydrate(smi):
    m = Chem.MolFromSmiles(smi)
    o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    c = o.GetNeighbors()[0]
    best = []
    for nb in c.GetNeighbors():
        if nb.GetSymbol() != 'C' or nb.GetTotalNumHs() == 0:
            continue
        rw = Chem.RWMol(m)
        rw.RemoveAtom(o.GetIdx())
        rw.GetBondBetweenAtoms(c.GetIdx(), nb.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
        mm = rw.GetMol()
        for a in mm.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
        Chem.SanitizeMol(mm)
        db = [b for b in mm.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE][0]
        nsub = sum(1 for at in (db.GetBeginAtom(), db.GetEndAtom()) for x in at.GetNeighbors()
                   if x.GetSymbol() == 'C' and x.GetIdx() not in (db.GetBeginAtomIdx(), db.GetEndAtomIdx()))
        best.append((nsub, Chem.MolToSmiles(mm)))
    best.sort(reverse=True)
    assert len(best) == 1 or best[0][0] > best[1][0] or best[0][1] == best[1][1]
    return best[0][1]

def tertiary(smi):
    m = Chem.MolFromSmiles(smi)
    o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    return o.GetNeighbors()[0].GetTotalNumHs() == 0

back, tert = [], []
for k, s in alk.items():
    assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C6H12'
    al = hydrate(s); ne = dehydrate(al)
    print(k, al, ne)
    if tertiary(al): tert.append(k)
    if Chem.CanonSmiles(ne) == Chem.CanonSmiles(s): back.append(k)
assert back == ['b', 'e'] and tert == ['a', 'b', 'd']
assert hydrate(alk['a']) == hydrate(alk['b'])
assert Chem.CanonSmiles(dehydrate(hydrate(alk['c']))) == Chem.CanonSmiles('CC=CC(C)C')
assert Chem.CanonSmiles(dehydrate(hydrate(alk['d']))) == Chem.CanonSmiles('CC(C)=C(C)C')
print('ORG-HYD-038 OK')
