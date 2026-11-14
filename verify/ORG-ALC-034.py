# ORG-ALC-034 検証：ケトンの α-H の重水素交換数とラセミ化（C7H14O 鎖式ケトンの総当たり）
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import product
def trees(n):
    for seq in product(range(n), repeat=n-2):
        deg = [1]*n
        for x in seq: deg[x] += 1
        edges = []
        for x in seq:
            leaf = min(i for i in range(n) if deg[i] == 1)
            edges.append((leaf, x)); deg[leaf] -= 1; deg[x] -= 1
        u, v = [i for i in range(n) if deg[i] == 1]; edges.append((u, v)); yield edges
# C7 の炭素骨格（木）に C=O を付ける（鎖の内部の C＝ケトン）
ket = set()
for edges in trees(7):
    em = Chem.RWMol()
    for i in range(7): em.AddAtom(Chem.Atom(6))
    for a, b in edges: em.AddBond(a, b, Chem.BondType.SINGLE)
    m = em.GetMol()
    if any(a.GetDegree() > 4 for a in m.GetAtoms()): continue
    for a in m.GetAtoms():
        if a.GetDegree() == 2:
            e2 = Chem.RWMol(m); o = e2.AddAtom(Chem.Atom(8)); e2.AddBond(a.GetIdx(), o, Chem.BondType.DOUBLE)
            p = e2.GetMol(); Chem.SanitizeMol(p); ket.add(Chem.MolToSmiles(p))
print('C7H14O 鎖式ケトン', len(ket))
rows = []
for s in sorted(ket):
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m) == 'C7H14O'
    co = [a for a in m.GetAtoms() if a.GetSymbol() == 'C' and any(n.GetSymbol() == 'O' for n in a.GetNeighbors())][0]
    alpha = [n for n in co.GetNeighbors() if n.GetSymbol() == 'C']
    nD = sum(a.GetTotalNumHs() for a in alpha)
    cc = [i for i, _ in Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)]
    rac = any(i in [a.GetIdx() for a in alpha if a.GetTotalNumHs() > 0] for i in cc)
    rows.append((s, bool(cc), rac, nD)); print(s, 'chiral' if cc else '', 'racemize' if rac else '', 'D+', nD)
P = [r for r in rows if r[1] and r[2] and r[3] == 3]
Q = [r for r in rows if r[1] and not r[2]]
R = [r for r in rows if r[1] and r[2] and r[3] == 4]
print('P', P, '\nQ', Q, '\nR', R)
assert [r[0] for r in P] == ['CCC(=O)C(C)CC']
assert [r[0] for r in Q] == ['CCC(C)CC(C)=O'] and Q[0][3] == 5
assert sorted(r[0] for r in R) == sorted(['CC(=O)C(C)C(C)C', 'CCCC(C)C(C)=O'])
# (1) 3-メチル-2-ブタノン
m = Chem.MolFromSmiles('CC(=O)C(C)C'); print('3-methyl-2-butanone alpha H =', 4)
# (4) 3-メチル-2-ペンタノンに対する 3,3-ジメチル（4級）では α に不斉 H なし → 光学活性の議論用
s5 = 'CCCC(C)(CC)C(C)=O'   # 3-エチル-3-メチル-2-ヘキサノン
m = Chem.MolFromSmiles(s5); print(CalcMolFormula(m))
cc = Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False); print('chiral', cc)
co = [a for a in m.GetAtoms() if a.GetSymbol()=='C' and any(n.GetSymbol()=='O' for n in a.GetNeighbors())][0]
alpha = [n for n in co.GetNeighbors() if n.GetSymbol()=='C']
assert len(cc) == 1 and cc[0][0] in [a.GetIdx() for a in alpha]
assert m.GetAtomWithIdx(cc[0][0]).GetTotalNumHs() == 0 and sum(a.GetTotalNumHs() for a in alpha) == 3
print('OK')
