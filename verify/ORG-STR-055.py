# ORG-STR-055 検証：C4H11N の異性体・アミド化・不斉炭素・C–H の種類
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import product
# 4 個の C と 1 個の N の鎖式飽和骨格をすべて列挙（木の SMILES を総当たり）
cands = set()
import itertools
atoms = ['C','C','C','C','N']
# 簡便法：知られている全骨格を明示列挙し，重複除去で 8 種を確認
smis = ['CCCCN','CCC(C)N','CC(C)CN','CC(C)(C)N','CCCNC','CC(C)NC','CCNCC','CCN(C)C',
        'NCCCC','CNCCC','CN(C)CC','CC(N)CC']
for s in smis:
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m) == 'C4H11N'
    cands.add(Chem.MolToSmiles(m))
print(len(cands), sorted(cands)); assert len(cands) == 8
# 総当たりでも確認：C4N の木（5 頂点）を全ラベル付き木から作る
from itertools import combinations
def trees(n):
    # Prüfer 列で全ラベル付き木
    for seq in product(range(n), repeat=n-2):
        deg = [1]*n
        for x in seq: deg[x] += 1
        edges = []; s = list(seq)
        for x in s:
            leaf = min(i for i in range(n) if deg[i] == 1)
            edges.append((leaf, x)); deg[leaf] -= 1; deg[x] -= 1
        u, v = [i for i in range(n) if deg[i] == 1]; edges.append((u, v)); yield edges
allc = set()
for edges in trees(5):
    for npos in range(5):
        em = Chem.RWMol()
        for i in range(5): em.AddAtom(Chem.Atom(7 if i == npos else 6))
        for a, b in edges: em.AddBond(a, b, Chem.BondType.SINGLE)
        m = em.GetMol()
        try:
            Chem.SanitizeMol(m)
        except Exception:
            continue
        if any(a.GetDegree() > 4 for a in m.GetAtoms()) or m.GetAtomWithIdx(npos).GetDegree() > 3: continue
        allc.add(Chem.MolToSmiles(m))
print('brute force', len(allc)); assert allc == cands
info = {}
for s in sorted(cands):
    m = Chem.MolFromSmiles(s); N = [a for a in m.GetAtoms() if a.GetSymbol() == 'N'][0]
    cls = N.GetDegree(); chiral = len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
    mh = Chem.AddHs(m); ranks = list(Chem.CanonicalRankAtoms(mh, breakTies=False))
    hk = {ranks[a.GetIdx()] for a in mh.GetAtoms() if a.GetSymbol() == 'H' and a.GetNeighbors()[0].GetSymbol() == 'C'}
    info[s] = (cls, chiral, len(hk)); print(s, '級', cls, '不斉', chiral, 'C-H 種類', len(hk))
ok = [s for s, (c, ch, k) in info.items() if c <= 2 and ch == 0 and k == 2]
print('D candidates', ok); assert ok == ['CCNCC']
assert sum(1 for v in info.values() if v[0] == 1) == 4 and sum(1 for v in info.values() if v[0] == 2) == 3
assert check_equation('(C2H5)2NH + (CH3CO)2O -> CH3CON(C2H5)2 + CH3COOH') is None
assert check_equation('(C2H5)2NH + HCl -> (C2H5)2NH2Cl') is None
print('OK')
