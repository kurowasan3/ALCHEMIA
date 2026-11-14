# ORG-STR-057 検証：C5H11NO2 の鎖式アミノ酸（-NH2 1個・-COOH 1個）の網羅と，COOH→H 置換後の不斉炭素
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D
RDLogger = __import__('rdkit').RDLogger; RDLogger.DisableLog('rdApp.*')
acids = ['CCCCC(=O)O', 'CCC(C)C(=O)O', 'CC(C)CC(=O)O', 'CC(C)(C)C(=O)O']   # C4H9-COOH の4骨格
res = set()
for a in acids:
    m = Chem.MolFromSmiles(a)
    for at in m.GetAtoms():
        if at.GetSymbol() == 'C' and at.GetTotalNumHs() > 0 and not any(n.GetSymbol() == 'O' for n in at.GetNeighbors()):
            rw = Chem.RWMol(m); i = rw.AddAtom(Chem.Atom(7)); rw.AddBond(at.GetIdx(), i, Chem.BondType.SINGLE)
            res.add(Chem.MolToSmiles(rw.GetMol()))
def nchiral(s):
    return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
rows = []
for s in sorted(res):
    m = Chem.MolFromSmiles(s); assert D.CalcMolFormula(m) == 'C5H11NO2'
    c = m.GetSubstructMatch(Chem.MolFromSmarts('C(=O)[OH]'))[0]
    n = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'N'][0]
    pos = {2: 'a', 3: 'b', 4: 'g', 5: 'd'}[len(Chem.GetShortestPath(m, c, n)) - 1]
    dec = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(Chem.ReplaceSubstructs(m, Chem.MolFromSmarts('C(=O)[OH]'), Chem.MolFromSmiles('[H]'))[0])))
    rows.append((s, pos, nchiral(s), dec, nchiral(dec)))
    print(rows[-1])
assert len(rows) == 12                                   # (1)
alpha = [r for r in rows if r[1] == 'a']; assert len(alpha) == 3   # (2)
assert all(r[2] == 1 for r in alpha)
keep = [r[0] for r in alpha if r[4] > 0]; assert keep == ['CCC(C)(N)C(=O)O']   # (3)
X = [r[0] for r in rows if r[1] == 'b' and r[2] > 0 and r[4] == 0]; assert X == ['CCC(CN)C(=O)O']  # (4)
print('OK')
