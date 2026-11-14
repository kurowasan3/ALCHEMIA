# Enumerate 5/6-membered hydroxy lactones C5H8O3 from C5H10O4 dihydroxy acids and apply clues
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D
from itertools import combinations
acids = ['OC(=O)CCCC', 'OC(=O)C(C)CC', 'OC(=O)CC(C)C', 'OC(=O)C(C)(C)C']
hyd = set()
for a in acids:
    m = Chem.MolFromSmiles(a)
    cs = [at.GetIdx() for at in m.GetAtoms() if at.GetSymbol() == 'C' and at.GetTotalNumHs() > 0
          and not any(n.GetSymbol() == 'O' for n in at.GetNeighbors())]
    for i, j in combinations(cs, 2):
        rw = Chem.RWMol(m)
        for k in (i, j):
            o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(k, o, Chem.BondType.SINGLE)
        x = rw.GetMol(); Chem.SanitizeMol(x); hyd.add(Chem.MolToSmiles(x))
assert len(hyd) == 17
def nst(m): return len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
acid = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
def lactones(h, sizes=(5, 6)):
    m = Chem.MolFromSmiles(h); c, _, oh = m.GetSubstructMatch(acid); out = []
    for a in m.GetAtoms():
        if a.GetSymbol() == 'O' and a.GetTotalNumHs() == 1 and a.GetIdx() != oh:
            ring = len(Chem.GetShortestPath(m, c, a.GetIdx()))
            if ring in sizes:
                rw = Chem.RWMol(m); rw.AddBond(c, a.GetIdx(), Chem.BondType.SINGLE); rw.RemoveAtom(oh)
                x = rw.GetMol(); Chem.SanitizeMol(x); out.append(Chem.MolToSmiles(x))
    return out
L = {}
for h in hyd:
    for s in lactones(h): L[s] = h
assert len(L) == 12 and all(D.CalcMolFormula(Chem.MolFromSmiles(s)) == 'C5H8O3' for s in L)
def oh_class(s):
    m = Chem.MolFromSmiles(s); i = m.GetSubstructMatch(Chem.MolFromSmarts('[CX4][OX2H1]'))[0]
    return sum(1 for n in m.GetAtomWithIdx(i).GetNeighbors() if n.GetSymbol() == 'C')
step = [s for s in L if oh_class(s) == 2 and nst(Chem.MolFromSmiles(s)) == 1]
print('after ア〜エ:', step); assert len(step) == 3
final = [s for s in step if any(x != s for x in lactones(L[s], sizes=(5,)))]
print('final:', final); assert final == [Chem.MolToSmiles(Chem.MolFromSmiles('O=C1CCC(O)CO1'))]
E = [x for x in lactones(L[final[0]], sizes=(5,))][0]; print('E', E)
B = 'O=C1CCC(=O)CO1'
print('stereo B', nst(Chem.MolFromSmiles(B)), 'E', nst(Chem.MolFromSmiles(E)))
assert nst(Chem.MolFromSmiles(B)) == 0 and nst(Chem.MolFromSmiles(E)) == 1
# combustion
n = 11.6/116; print('CO2 mg', n*5*44, 'H2O mg', n*4*18)
assert abs(n*5*44 - 22.0) < 1e-9 and abs(n*4*18 - 7.20) < 1e-9
print('OK')
