# C5H10O carbonyl compounds: reduction -> alcohol; dehydration -> alkenes (all beta-H, no rearrangement)
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
cands = {'pentanal':'CCCCC=O','2-methylbutanal':'CCC(C)C=O','3-methylbutanal':'CC(C)CC=O',
 '2,2-dimethylpropanal':'CC(C)(C)C=O','2-pentanone':'CCCC(C)=O','3-pentanone':'CCC(=O)CC','3-methyl-2-butanone':'CC(C)C(C)=O'}
def nchiral(m): return len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
rxn_red = None
def reduce(smi):
    m = Chem.MolFromSmiles(smi); rw = Chem.RWMol(m)
    for b in rw.GetBonds():
        if b.GetBondType()==Chem.BondType.DOUBLE: b.SetBondType(Chem.BondType.SINGLE)
    s = Chem.MolToSmiles(rw); return Chem.MolToSmiles(Chem.MolFromSmiles(s))
def dehydrate(smi):
    m = Chem.MolFromSmiles(smi); out=set()
    o = [a for a in m.GetAtoms() if a.GetSymbol()=='O'][0]; ca = o.GetNeighbors()[0]
    for nb in ca.GetNeighbors():
        if nb.GetSymbol()=='C' and nb.GetTotalNumHs()>0:
            rw = Chem.RWMol(m); rw.GetBondBetweenAtoms(ca.GetIdx(), nb.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
            rw.RemoveAtom(o.GetIdx()); x = rw.GetMol(); Chem.SanitizeMol(x); out.add(Chem.MolToSmiles(x))
    return out
def cistrans(smi):
    m = Chem.MolFromSmiles(smi); return len(list(EnumerateStereoisomers(m)))>1
res = {}
for n,s in cands.items():
    m = Chem.MolFromSmiles(s); assert D.CalcMolFormula(m)=='C5H10O'
    al = reduce(s); alk = dehydrate(al)
    io = m.HasSubstructMatch(Chem.MolFromSmarts('[CH3]C(=O)[#6,#1]')) or m.HasSubstructMatch(Chem.MolFromSmarts('[CH3]C=O'))
    ald = m.HasSubstructMatch(Chem.MolFromSmarts('[CX3H1](=O)'))
    res[n] = (nchiral(m), nchiral(Chem.MolFromSmiles(al)), sorted(alk), [cistrans(a) for a in sorted(alk)], io, ald)
    print(n, res[n])
# (1) ketone whose alcohol has no stereocenter
assert [n for n,r in res.items() if not r[5] and r[1]==0] == ['3-pentanone']
# (2) carbonyl whose alcohol cannot dehydrate
assert [n for n,r in res.items() if len(r[2])==0] == ['2,2-dimethylpropanal']
# (3) chiral itself and alcohol chiral count same
assert [n for n,r in res.items() if r[0]==1] == ['2-methylbutanal']
# (4) iodoform-positive ketone, alkenes all without cis-trans
print([n for n,r in res.items() if r[4] and not r[5] and not any(r[3])])
assert [n for n,r in res.items() if r[4] and not r[5] and not any(r[3])] == ['3-methyl-2-butanone']
# (5) compounds whose reduction increases stereocenters
print([n for n,r in res.items() if r[1]>r[0]])
print('OK')
