# ORG-ARO-081 検証：C10H14 芳香族異性体の網羅と構造決定
from rdkit import Chem
from itertools import combinations, product
from collections import Counter
alk = {1: ['C'], 2: ['CC'], 3: ['CCC', 'C(C)C'], 4: ['CCCC', 'C(C)CC', 'CC(C)C', 'C(C)(C)C']}
def parts(n, k):
    if k == 1: yield (n,); return
    for i in range(1, n):
        for p in parts(n - i, k - 1): yield (i,) + p
iso = set()
for k in range(1, 5):
    for pos in combinations(range(6), k):
        for sizes in parts(4, k):
            for grp in product(*[alk[s] for s in sizes]):
                at = []
                for i in range(6):
                    sub = '(' + grp[pos.index(i)] + ')' if i in pos else ''
                    at.append(('c1' if i == 0 else 'c') + sub)
                smi = ''.join(at[:5]) + 'c' + ('(' + grp[pos.index(5)] + ')' if 5 in pos else '') + '1'
                m = Chem.MolFromSmiles(smi)
                assert m is not None
                iso.add(Chem.MolToSmiles(m))
assert len(iso) == 22
def ncls(m):
    r = list(Chem.CanonicalRankAtoms(m, breakTies=False)); return len({r[a.GetIdx()] for a in m.GetAtoms()})
def oxid(m):
    # 環に結合した炭素ごとに，その炭素に H があれば COOH
    n = 0; pos = []
    for a in m.GetAtoms():
        if not a.GetIsAromatic() and any(nb.GetIsAromatic() for nb in a.GetNeighbors()):
            if a.GetTotalNumHs() > 0:
                n += 1; pos.append([nb.GetIdx() for nb in a.GetNeighbors() if nb.GetIsAromatic()][0])
    return n, pos
def ortho(m, pos):
    return any(m.GetBondBetweenAtoms(i, j) for i, j in combinations(pos, 2))
def para(m, pos):
    return len(pos) == 2 and len(Chem.GetShortestPath(m, pos[0], pos[1])) == 4
D = {s: Chem.MolFromSmiles(s) for s in iso}
print(sorted(Counter(ncls(m) for m in D.values()).items()))
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
P = [s for s, m in D.items() if ncls(m) == 3]; assert P == [can('Cc1cc(C)c(C)cc1C')]
Qc = [s for s, m in D.items() if ncls(m) == 4]; assert Qc == [can('CCc1ccc(CC)cc1')]
R = [s for s, m in D.items() if ncls(m) == 5 and oxid(m)[0] == 2]; assert R == [can('CCc1ccccc1CC')]
S = [s for s, m in D.items() if ncls(m) == 6 and oxid(m)[0] == 0]; assert S == [can('CC(C)(C)c1ccccc1')]
T = [s for s, m in D.items() if ncls(m) == 7 and oxid(m)[0] == 3 and ortho(m, oxid(m)[1])]; assert T == [can('CCc1c(C)cccc1C')]
T2 = [s for s, m in D.items() if ncls(m) == 7 and oxid(m)[0] == 3 and not ortho(m, oxid(m)[1])]; assert T2 == [can('CCc1cc(C)cc(C)c1')]
U = [s for s, m in D.items() if ncls(m) == 7 and oxid(m)[0] == 2 and para(m, oxid(m)[1])]; assert U == [can('Cc1ccc(C(C)C)cc1')]
mono = [s for s, m in D.items() if len([a for a in m.GetAtoms() if a.GetIsAromatic() and a.GetTotalNumHs() == 0]) == 1]
assert len(mono) == 4 and [s for s in mono if oxid(D[s])[0] == 0] == [can('CC(C)(C)c1ccccc1')]
ten = [s for s, m in D.items() if ncls(m) == 10]; print(ten); assert len(ten) == 6
print('OK')
