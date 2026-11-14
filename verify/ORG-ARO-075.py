# ORG-ARO-075 検証：C8H10 の非等価な炭素と一塩素置換体
from rdkit import Chem
def classes(smi, pred):
    m = Chem.MolFromSmiles(smi); r = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    return len({r[a.GetIdx()] for a in m.GetAtoms() if pred(a)})
def ring_cl(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    for a in m.GetAtoms():
        if a.GetIsAromatic() and a.GetTotalNumHs() == 1:
            rw = Chem.RWMol(Chem.AddHs(m))
            h = [n.GetIdx() for n in rw.GetAtomWithIdx(a.GetIdx()).GetNeighbors() if n.GetSymbol() == 'H'][0]
            rw.GetAtomWithIdx(h).SetAtomicNum(17)
            out.add(Chem.MolToSmiles(Chem.RemoveHs(rw)))
    return len(out)
D = {'CCc1ccccc1': (6, 3), 'Cc1ccccc1C': (4, 2), 'Cc1cccc(C)c1': (5, 3), 'Cc1ccc(C)cc1': (3, 1)}
for s, (nc, ncl) in D.items():
    a = classes(s, lambda x: x.GetSymbol() == 'C'); b = ring_cl(s); print(s, a, b)
    assert (a, b) == (nc, ncl)
assert classes('OC(=O)c1ccccc1', lambda x: x.GetSymbol() == 'C') == 5
print('OK')
