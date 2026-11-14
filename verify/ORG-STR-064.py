# ORG-STR-064 検証：ジベンゾ-p-ジオキシン（P）・ジベンゾフラン（Q）の塩素置換体の数（RDKit の正規化 SMILES で重複除去）
from itertools import combinations
from rdkit import Chem
def subst(smi, k):
    m = Chem.MolFromSmiles(smi)
    Hs = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'C' and a.GetTotalNumHs() > 0]
    out = {}
    for c in combinations(Hs, k):
        rw = Chem.RWMol(m)
        for i in c:
            j = rw.AddAtom(Chem.Atom(17)); rw.AddBond(i, j, Chem.BondType.SINGLE)
        mm = rw.GetMol(); Chem.SanitizeMol(mm)
        out.setdefault(Chem.MolToSmiles(mm), []).append(c)
    return out
P = 'c1ccc2c(c1)Oc1ccccc1O2'      # 原子番号 → 位置番号：2→1, 1→2, 0→3, 5→4, 8→6, 9→7, 10→8, 11→9
Q = 'c1ccc2c(c1)oc1ccccc12'
pos = {2: 1, 1: 2, 0: 3, 5: 4, 8: 6, 9: 7, 10: 8, 11: 9}
n1P, n1Q = len(subst(P, 1)), len(subst(Q, 1))
n2P, n2Q = len(subst(P, 2)), len(subst(Q, 2))
print('mono P,Q', n1P, n1Q, 'di P,Q', n2P, n2Q)
assert (n1P, n1Q, n2P, n2Q) == (2, 4, 10, 16)
ringA = {1, 2, 3, 4}
same = diff = 0
names = {}
for smi, combos in subst(P, 2).items():
    c = sorted(pos[i] for i in combos[0])
    names[smi] = sorted(tuple(sorted(pos[i] for i in cc)) for cc in combos)
    if (c[0] in ringA) == (c[1] in ringA): same += 1
    else: diff += 1
print('same ring', same, 'different rings', diff)
assert (same, diff) == (4, 6)
# 1,6 と 1,9 は別物，2,7 と 2,8 は別物
def smi_of(a, b):
    inv = {v: k for k, v in pos.items()}
    rw = Chem.RWMol(Chem.MolFromSmiles(P))
    for p in (a, b):
        j = rw.AddAtom(Chem.Atom(17)); rw.AddBond(inv[p], j, Chem.BondType.SINGLE)
    m = rw.GetMol(); Chem.SanitizeMol(m); return Chem.MolToSmiles(m)
assert smi_of(1, 6) != smi_of(1, 9) and smi_of(2, 7) != smi_of(2, 8)
assert smi_of(1, 6) == smi_of(4, 9) and smi_of(1, 9) == smi_of(4, 6)
# 同じ環に Cl 2 個の 4 種について，三塩素化物の数
tri = {}
for pair in [(1, 2), (1, 3), (1, 4), (2, 3)]:
    tri[pair] = len(subst(smi_of(*pair), 1))
print('trichloro products', tri)
assert tri == {(1, 2): 6, (1, 3): 6, (1, 4): 3, (2, 3): 3}
X = smi_of(2, 3)
n4 = len(subst(X, 2))
print('tetrachloro from 2,3-dichloro', n4, 'total tetra', len(subst(P, 4)))
assert n4 == 9
assert smi_of(2, 3) != smi_of(1, 4)
# 2,3,7,8-四塩素化物が含まれる
inv = {v: k for k, v in pos.items()}
rw = Chem.RWMol(Chem.MolFromSmiles(P))
for p in (2, 3, 7, 8):
    j = rw.AddAtom(Chem.Atom(17)); rw.AddBond(inv[p], j, Chem.BondType.SINGLE)
m = rw.GetMol(); Chem.SanitizeMol(m)
assert Chem.MolToSmiles(m) in subst(X, 2)
print('OK')
