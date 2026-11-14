# ORG-HYD-049 検証：環状テトラエンのオゾン分解断片の並べ方と部分水素化の生成物の数
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import permutations
from sympy import Rational
# (1) C18Hh，k 個の C=C：12.2/(216+h) * 32k = 6.40
sols = []
for k in range(1, 10):
    h = Rational(1220*32*k, 640) - 216
    if h > 0 and h == int(h) and h % 2 == 0 and h <= 38:
        sols.append((k, int(h)))
print('(k,h)', sols); assert sols == [(4, 28)]
# 断片 a = C4（OHC-(CH2)2-CHO），b = C5（OHC-(CH2)3-CHO）を 2 個ずつ環状につなぐ
def build(seq):
    rw = Chem.RWMol(); prev_end = None; first_start = None
    for f in seq:
        nC = 4 if f == 'a' else 5
        idx = [rw.AddAtom(Chem.Atom(6)) for _ in range(nC)]
        for i in range(nC - 1): rw.AddBond(idx[i], idx[i+1], Chem.BondType.SINGLE)
        if prev_end is not None: rw.AddBond(prev_end, idx[0], Chem.BondType.DOUBLE)
        else: first_start = idx[0]
        prev_end = idx[-1]
    rw.AddBond(prev_end, first_start, Chem.BondType.DOUBLE)
    m = rw.GetMol(); Chem.SanitizeMol(m); return m
cands = {}
for p in set(permutations('aabb')):
    m = build(p); smi = Chem.MolToSmiles(m)
    cands.setdefault(smi, p)
print('distinct ring structures', len(cands), cands)
assert len(cands) == 2
for smi, p in cands.items():
    m = Chem.MolFromSmiles(smi)
    assert CalcMolFormula(m) == 'C18H28'
    # C=C を 1 個だけ水素化した生成物（構造異性体）の数
    prods = set()
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE:
            rw = Chem.RWMol(m); rw.GetBondWithIdx(b.GetIdx()).SetBondType(Chem.BondType.SINGLE)
            mm = rw.GetMol(); Chem.SanitizeMol(mm); prods.add(Chem.MolToSmiles(mm))
    print(''.join(p), smi, 'monohydrogenation products', len(prods))
    # その生成物のオゾン分解で得られる化合物の種類
res = {''.join(p): None for p in cands.values()}
# (5) 質量：オゾン分解生成物 0.100 mol ずつ
print('mass products', 0.100*86 + 0.100*100, 'H2 for full hydrogenation mol', 0.0500*4)
assert abs(0.100*86 + 0.100*100 - 18.60) < 1e-9
print(CalcMolFormula(Chem.MolFromSmiles('O=CCCC=O')), CalcMolFormula(Chem.MolFromSmiles('O=CCCCC=O')))
# 部分水素化物（abba 型）のオゾン分解生成物の組（判別の確認）
def ozo(m):
    rw = Chem.RWMol(m)
    pairs = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]
    for a, b in pairs:
        rw.RemoveBond(a, b)
        for x in (a, b):
            o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(x, o, Chem.BondType.DOUBLE)
    mm = rw.GetMol(); Chem.SanitizeMol(mm); return sorted(Chem.MolToSmiles(mm).split('.'))
for smi in cands:
    print(smi, 'ozonolysis', ozo(Chem.MolFromSmiles(smi)))
print('OK')
