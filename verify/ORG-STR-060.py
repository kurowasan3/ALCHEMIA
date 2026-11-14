from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import numpy as np
RDLogger.DisableLog('rdApp.*')
def sym_classes(smi):
    """3D 構造の局所的な向き（4 配位原子の符号付き体積）を保つ（回転）か，すべて反転する（鏡映）グラフ自己同型だけを分子の対称操作とみなし，炭素の軌道（等価な組）を数える"""
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=7)
    X = m.GetConformer().GetPositions()
    four = [a.GetIdx() for a in m.GetAtoms() if a.GetDegree() == 4]
    nb = {i: [n.GetIdx() for n in m.GetAtomWithIdx(i).GetNeighbors()] for i in four}
    def vol(i, order):
        a, b, c, d = (X[k] - X[i] for k in order); return np.sign(np.dot(a - d, np.cross(b - d, c - d)))
    sgn = {i: vol(i, nb[i]) for i in four}
    good = []; improper = False
    for s in m.GetSubstructMatches(m, uniquify=False, useChirality=False, maxMatches=2000000):
        rel = set()
        for i in four:
            j = s[i]; inv = {s[k]: k for k in nb[i]}
            # j の近接原子を，i の近接原子の像として並べた順序での符号
            img = [s[k] for k in nb[i]]
            rel.add(vol(j, img) == sgn[i])
            if len(rel) > 1: break
        if len(rel) > 1: continue
        good.append(s)
        if rel == {False}: improper = True
    orb = {frozenset(s[i] for s in good) for i, a in enumerate(m.GetAtoms()) if a.GetSymbol() == 'C'}
    return len(orb), improper   # improper=True なら鏡映対称をもつ（キラルでない）
def stereoisomers(smi):
    m = Chem.MolFromSmiles(smi)
    return sorted({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m, options=StereoEnumerationOptions(unique=True, onlyUnassigned=True))})
# ---- ORG-STR-060 検証 ----
from itertools import combinations
def grow(smis, z):
    out = set()
    for s in smis:
        m = Chem.MolFromSmiles(s)
        for a in m.GetAtoms():
            rw = Chem.RWMol(m); i = rw.AddAtom(Chem.Atom(z)); rw.AddBond(a.GetIdx(), i, Chem.BondType.SINGLE)
            mm = rw.GetMol()
            try: Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
            except Exception: pass
    return out
lv = {'C', 'O'}
for k in range(5): lv = grow(lv, 6) | grow(lv, 8)
trees = {t for t in lv if Chem.MolFromSmiles(t).GetNumAtoms() == 6 and sum(a.GetSymbol() == 'O' for a in Chem.MolFromSmiles(t).GetAtoms()) == 1}
iso = set()
for s in trees:
    m = Chem.MolFromSmiles(s)
    for i, j in combinations(range(m.GetNumAtoms()), 2):
        rw = Chem.RWMol(m); b = rw.GetBondBetweenAtoms(i, j)
        if b is None: rw.AddBond(i, j, Chem.BondType.SINGLE)
        else: b.SetBondType(Chem.BondType.DOUBLE)
        mm = rw.GetMol()
        try: Chem.SanitizeMol(mm); iso.add(Chem.MolToSmiles(mm))
        except Exception: pass
print('C5H10O の構造異性体', len(iso))
cyc = sorted(s for s in iso if '=' not in s)
print('環をもち二重結合のないもの', len(cyc))
def nC_flat(s):
    m = Chem.MolFromSmiles(s); r = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    return len({r[a.GetIdx()] for a in m.GetAtoms() if a.GetSymbol() == 'C'})
c3 = [s for s in cyc if nC_flat(s) == 3]; print('(1) 構造式上で3種類:', c3)
assert sorted(c3) == sorted(['C1CCOCC1', 'CC1(C)COC1', 'CC1C(C)C1O', 'CC1CC(C)O1', 'OC1CCCC1'])
# 立体異性体ごと（他の 28 種類も含め，立体で 3 種類になるものがないか確認）
three = []
for s in cyc:
    for i in stereoisomers(s):
        n, achiral = sym_classes(i)
        if n == 3: three.append((s, i, achiral))
for t in three: print('  3種類:', t)
assert len(three) == 8
# (2) 2,3-ジメチルシクロプロパノール
st = stereoisomers('CC1C(C)C1O'); res = [sym_classes(i) for i in st]; print('(2)', list(zip(st, res)))
assert len(st) == 4 and sorted(r[0] for r in res) == [3, 3, 5, 5]
assert all(r[1] for r in res if r[0] == 3) and not any(r[1] for r in res if r[0] == 5)
hs = Chem.FindMolChiralCenters(Chem.MolFromSmiles('CC1C(C)C1O'), includeUnassigned=True, useLegacyImplementation=True)
print('  高校の定義の不斉炭素（4つの異なる原子団）:', hs)
# (3) 2,4-ジメチルオキセタン
st = stereoisomers('CC1CC(C)O1'); res = [sym_classes(i) for i in st]; print('(3)', list(zip(st, res)))
assert len(st) == 3 and all(r[0] == 3 for r in res) and sum(r[1] for r in res) == 1
# (4) 酸化生成物の不斉炭素
for k in ['O=C1CCCC1', 'CC1C(C)C1=O']:
    print('(4)', k, Chem.FindMolChiralCenters(Chem.MolFromSmiles(k), includeUnassigned=True, useLegacyImplementation=True))
print('OK')
