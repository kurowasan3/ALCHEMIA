ITEM = '/home/claude/chem2/outh/items/ORG-STR-020/item.json'
# 検証スクリプト（グループ g08）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return eqs

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def check_parts():
    n_items = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, n_items + 1)], (n_items, nums)
    print('小問の数OK:', n_items)


check_parts()
import itertools
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
RDLogger.DisableLog('rdApp.*')
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
BT = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
iso = set(); pairs = list(itertools.combinations(range(4), 2))
for orders in itertools.product(range(4), repeat=6):
    if sum(orders) != 5: continue
    rw = Chem.RWMol()
    for _ in range(4): rw.AddAtom(Chem.Atom(6))
    for (a, b), o in zip(pairs, orders):
        if o: rw.AddBond(a, b, BT[o])
    m = rw.GetMol()
    try: Chem.SanitizeMol(m)
    except Exception: continue
    if len(Chem.GetMolFrags(m)) == 1: iso.add(Chem.MolToSmiles(m))
assert len(iso) == 9, iso
def rings(s): return Chem.MolFromSmiles(s).GetRingInfo().NumRings()
def pi(s): return sum({1:0, 2:1, 3:2}[int(b.GetBondTypeAsDouble())] for b in Chem.MolFromSmiles(s).GetBonds())
cls = {}
for s in iso: cls.setdefault((rings(s), pi(s)), []).append(s)
print({k: len(v) for k, v in cls.items()}); assert {k: len(v) for k, v in cls.items()} == {(0, 2): 4, (1, 1): 4, (2, 0): 1}
def hyd(s):
    m = Chem.RWMol(Chem.MolFromSmiles(s))
    for b in m.GetBonds(): b.SetBondType(Chem.BondType.SINGLE)
    for a in m.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    mm = m.GetMol(); Chem.SanitizeMol(mm); return Chem.MolToSmiles(mm)
def nst(s): return len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(s), options=StereoEnumerationOptions(unique=True))})
P = [s for s in iso if pi(s) == 1 and hyd(s) == can('C1CCC1')]; assert P == [can('C1=CCC1')]
R = [s for s in iso if pi(s) == 0]; assert R == [can('C1C2CC12')]
adduct = {can('CC1=CC1'): 'CC1(Br)CC1Br', can('CC1C=C1'): 'CC1C(Br)C1Br', can('C=C1CC1'): 'BrCC1(Br)CC1'}
Qc = [s for s in iso if pi(s) == 1 and hyd(s) == can('CC1CC1')]; assert sorted(Qc) == sorted(adduct)
cnt = {s: nst(adduct[s]) for s in Qc}; print(cnt)
Q = [s for s in Qc if cnt[s] == 1]; assert Q == [can('C=C1CC1')] and cnt[can('CC1=CC1')] == 4 and cnt[can('CC1C=C1')] == 4
SC = [s for s in iso if pi(s) == 2 and hyd(s) == can('CCCC')]; assert len(SC) == 4
term = Chem.MolFromSmarts('[CH]#C'); S = [s for s in SC if Chem.MolFromSmiles(s).HasSubstructMatch(term)]; assert S == [can('C#CCC')]
def linear4(s):
    m = Chem.MolFromSmiles(s)
    sp = [a.GetHybridization() == Chem.HybridizationType.SP for a in m.GetAtoms()]
    # 4 個の C が一直線：内側の 2 個がともに sp
    inner = [a.GetIdx() for a in m.GetAtoms() if a.GetDegree() == 2]
    return len(inner) == 2 and all(sp[i] for i in inner)
T = [s for s in SC if linear4(s)]; assert T == [can('CC#CC')], T
# (4) 自身の立体異性体
counts = {s: nst(s) for s in iso}; print(counts)
for s, n in counts.items():
    if s == can('C1C2CC12'):
        forms = {Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(s))}
        ok = []
        for f in forms:
            good = False
            for seed in range(1, 6):
                m = Chem.AddHs(Chem.MolFromSmiles(f)); p = AllChem.ETKDGv3(); p.useRandomCoords = True; p.randomSeed = seed; p.maxIterations = 2000
                if AllChem.EmbedMolecule(m, p) == 0:
                    good = True; break
            ok.append(good)
        assert sum(ok) == 1, ok
        real = [f for f, k in zip(forms, ok) if k][0]
        mir = real.replace('@@', '!').replace('@', '@@').replace('!', '@')
        assert can(mir) == can(real)  # 存在できる形はアキラル
    else:
        assert n == 1, s

print('全検証OK:', ITEM)
