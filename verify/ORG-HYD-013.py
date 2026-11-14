# ---- small isomer enumerator (RDKit) ----
from rdkit import Chem, RDLogger
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
RDLogger.DisableLog('rdApp.*')
VAL = {'C': 4, 'O': 2, 'N': 3, 'Cl': 1, 'Br': 1}

def _canon(rw):
    m = rw.GetMol()
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    return Chem.MolToSmiles(m)

def enumerate_isomers(heavy, nH):
    """heavy: dict like {'C':4,'O':1}; returns set of canonical SMILES (constitutional isomers)."""
    from collections import Counter
    target = Counter(heavy)
    total = sum(target.values())
    first = [e for e in target][0]
    start = Chem.RWMol(); start.AddAtom(Chem.Atom(first))
    frontier = {_canon(start)}
    seen = set(frontier)
    results = set()
    while frontier:
        new = set()
        for smi in frontier:
            m = Chem.RWMol(Chem.MolFromSmiles(smi))
            Chem.Kekulize(m, clearAromaticFlags=True)
            have = Counter(a.GetSymbol() for a in m.GetAtoms())
            def free(i):
                a = m.GetAtomWithIdx(i)
                used = sum(b.GetBondTypeAsDouble() for b in a.GetBonds())
                return VAL[a.GetSymbol()] - used
            n = m.GetNumAtoms()
            cand = []
            # add atom
            if sum(have.values()) < total:
                for e in target:
                    if have[e] >= target[e]:
                        continue
                    for i in range(n):
                        for order in (1, 2, 3):
                            if free(i) >= order and VAL[e] >= order:
                                r = Chem.RWMol(m)
                                j = r.AddAtom(Chem.Atom(e))
                                r.AddBond(i, j, {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[order])
                                cand.append(r)
            # raise order / add ring bond
            for i in range(n):
                for j in range(i + 1, n):
                    if free(i) >= 1 and free(j) >= 1:
                        b = m.GetBondBetweenAtoms(i, j)
                        r = Chem.RWMol(m)
                        if b is None:
                            r.AddBond(i, j, Chem.BondType.SINGLE)
                        else:
                            o = b.GetBondTypeAsDouble()
                            if o >= 3:
                                continue
                            r.GetBondBetweenAtoms(i, j).SetBondType(Chem.BondType.DOUBLE if o == 1 else Chem.BondType.TRIPLE)
                        cand.append(r)
            for r in cand:
                for a in r.GetAtoms():
                    a.SetNoImplicit(False)
                c = _canon(r)
                if c is None or c in seen:
                    continue
                mm = Chem.AddHs(Chem.MolFromSmiles(c))
                cnt = Counter(a.GetSymbol() for a in mm.GetAtoms())
                h = cnt.pop('H', 0)
                rem = target - cnt
                if h + sum((VAL[e] - 2) * k for e, k in rem.items()) < nH:
                    continue
                seen.add(c); new.add(c)
                if cnt == target and h == nH:
                    results.add(c)
        frontier = new
    return results

def stereo_isomers(smi):
    m = Chem.MolFromSmiles(smi)
    opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True)
    return sorted({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m, options=opts)})

def has(smi, smarts):
    return Chem.MolFromSmiles(smi).HasSubstructMatch(Chem.MolFromSmarts(smarts))

def n_chiral(smi):
    m = Chem.MolFromSmiles(smi)
    return len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))

def canon(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))
# ---- end enumerator ----

# ORG-HYD-013 検証：マルコフニコフ則による水の付加の主生成物・副生成物
ring = Chem.MolFromSmarts('[R]')
alk = [s for s in enumerate_isomers({'C': 5}, 10) if not Chem.MolFromSmiles(s).HasSubstructMatch(ring)]
def add(smi, X=8):
    m = Chem.MolFromSmiles(smi); b = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE][0]
    a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
    out = {}
    for tgt, other in ((a1, a2), (a2, a1)):
        rw = Chem.RWMol(m); rw.GetBondBetweenAtoms(a1.GetIdx(), a2.GetIdx()).SetBondType(Chem.BondType.SINGLE)
        j = rw.AddAtom(Chem.Atom(X)); rw.AddBond(tgt.GetIdx(), j, Chem.BondType.SINGLE)
        for a in rw.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
        mm = rw.GetMol(); Chem.SanitizeMol(mm)
        # マルコフニコフ：H は H の多い炭素に → X は H の少ない炭素に
        out[canon(Chem.MolToSmiles(mm))] = tgt.GetTotalNumHs() - other.GetTotalNumHs()
    vals = sorted(set(out.values()))
    if len(out) == 1: return list(out), []
    if vals[0] == vals[-1]: return sorted(out), []        # 同程度
    major = [k for k, v in out.items() if v == vals[0]]; minor = [k for k, v in out.items() if v == vals[-1]]
    return major, minor
def branched(s):
    m = Chem.MolFromSmiles(s); return any(sum(n.GetSymbol() == 'C' for n in a.GetNeighbors()) >= 3 for a in m.GetAtoms())
def cls(alc):
    m = Chem.MolFromSmiles(alc); c = m.GetSubstructMatch(Chem.MolFromSmarts('[CX4][OH]'))[0]
    return {2: 1, 1: 2, 0: 3}[m.GetAtomWithIdx(c).GetTotalNumHs()] if m.GetAtomWithIdx(c).GetTotalNumHs() < 3 else 0
br = [s for s in alk if branched(s)]
info = {s: add(s) for s in br}
for k, v in info.items(): print(k, v)
tert = [s for s, (M, m) in info.items() if len(M) == 1 and cls(M[0]) == 3]
D = {info[s][0][0] for s in tert}; assert len(D) == 1 and D == {canon('CCC(C)(C)O')}
A = [s for s in tert if cls(info[s][1][0]) == 2]; B = [s for s in tert if cls(info[s][1][0]) == 1]
C = [s for s in br if s not in tert]
assert A == [canon('CC=C(C)C')] and B == [canon('C=C(C)CC')] and C == [canon('C=CC(C)C')]
E, F = info[C[0]][0][0], info[C[0]][1][0]
assert E == canon('CC(O)C(C)C') and F == canon('CC(C)CCO')
# E の酸化 → 3-メチル-2-ブタノン（CH3CO- あり），F の酸化 → アルデヒド
assert has('CC(=O)C(C)C', '[CH3]C(=O)[#6]')
# 2-ペンテン：両端の炭素の H はともに1個 → 2-ペンタノールと 3-ペンタノールが同程度
M, m = add('CC=CCC'); print('2-pentene ->', M, m); assert sorted(M) == sorted([canon('CCCC(C)O'), canon('CCC(O)CC')]) and m == []
# HBr の付加（A）
assert add(A[0], 35)[0] == [canon('CCC(C)(C)Br')]
print('ORG-HYD-013 OK')
