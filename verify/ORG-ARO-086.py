# ORG-ARO-086 検証：一置換ベンゼン C9H11NO のアミドの網羅と加水分解生成物による分類
import itertools
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, FindMolChiralCenters
def prufer_trees(n):
    for seq in itertools.product(range(n), repeat=n - 2):
        deg = [1] * n
        for s in seq: deg[s] += 1
        edges = []
        for s in seq:
            for i in range(n):
                if deg[i] == 1:
                    edges.append((i, s)); deg[i] -= 1; deg[s] -= 1; break
        edges.append(tuple(i for i in range(n) if deg[i] == 1))
        yield edges
# 頂点 0=フェニル，1-3=C，4=N。側鎖は鎖式で C=O 以外に不飽和なし，C=O の炭素は N と結合（アミド）
res = {}
for edges in prufer_trees(5):
    deg = [0] * 5
    for a, b in edges: deg[a] += 1; deg[b] += 1
    if deg[0] != 1 or deg[4] > 3: continue
    for co in (1, 2, 3):
        if not ((co, 4) in edges or (4, co) in edges) or deg[co] > 2: continue
        rw = Chem.RWMol(Chem.MolFromSmiles('c1ccccc1'))
        idx = {0: 0}
        for k, sym in [(1, 'C'), (2, 'C'), (3, 'C'), (4, 'N')]:
            idx[k] = rw.AddAtom(Chem.Atom(sym))
        for a, b in edges: rw.AddBond(idx[a], idx[b], Chem.BondType.SINGLE)
        o = rw.AddAtom(Chem.Atom('O')); rw.AddBond(idx[co], o, Chem.BondType.DOUBLE)
        m = rw.GetMol(); Chem.SanitizeMol(m)
        res[Chem.MolToSmiles(m)] = m
assert all(rdMolDescriptors.CalcMolFormula(m) == 'C9H11NO' for m in res.values())
print('total', len(res)); assert len(res) == 12
rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])[N:3]>>[C:1](=[O:2])O.[N:3]')
benz = Chem.MolFromSmarts('c1ccccc1')
cnt = dict(chiral=0, arom_acid=0, formic=0, tert=0)
amines = {}
for s, m in sorted(res.items()):
    if FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False): cnt['chiral'] += 1
    acid, amine = rxn.RunReactants((m,))[0]
    Chem.SanitizeMol(acid); Chem.SanitizeMol(amine)
    a, n = Chem.MolToSmiles(acid), Chem.MolToSmiles(amine)
    if acid.HasSubstructMatch(benz): cnt['arom_acid'] += 1
    if a == 'O=CO': cnt['formic'] += 1
    namide = m.GetAtomWithIdx([x.GetIdx() for x in m.GetAtoms() if x.GetSymbol() == 'N'][0])
    if namide.GetTotalNumHs() == 0: cnt['tert'] += 1
    amines[s] = (a, n)
    print(f'{s:24s} acid {a:18s} amine {n}')
print(cnt)
assert cnt == dict(chiral=2, arom_acid=5, formic=4, tert=4)
# 構造決定：酸がギ酸でない，アミンがベンゼン環をもつ第二級アミン
cand = []
for s, (a, n) in amines.items():
    am = Chem.MolFromSmiles(n)
    N = [x for x in am.GetAtoms() if x.GetSymbol() == 'N'][0]
    if a != 'O=CO' and am.HasSubstructMatch(benz) and N.GetDegree() == 2:
        cand.append(s)
print('Y', cand)
assert cand == [Chem.CanonSmiles('CC(=O)N(C)c1ccccc1')]
# アニリン（ジアゾ化できる第一級芳香族アミン）を生じるもの
anil = [s for s, (a, n) in amines.items() if n == 'Nc1ccccc1']
assert anil == [Chem.CanonSmiles('CCC(=O)Nc1ccccc1')]
print('OK')
