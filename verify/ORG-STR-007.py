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

# ORG-STR-007 検証：異性体の数（RDKit で全構造を列挙）
ring = Chem.MolFromSmarts('[R]')
c5h10 = enumerate_isomers({'C': 5}, 10)
chain = [s for s in c5h10 if not Chem.MolFromSmiles(s).HasSubstructMatch(ring)]
chain_st = [x for s in chain for x in stereo_isomers(s)]
print('C5H10 chain:', len(chain), 'with stereo:', len(chain_st), chain_st)
assert len(chain_st) == 6
# 水素付加後のアルカン
from collections import defaultdict
grp = defaultdict(list)
for s in chain_st:
    m = Chem.MolFromSmiles(s)
    h = Chem.RWMol(m)
    for b in h.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetStereo(Chem.BondStereo.STEREONONE)
    grp[canon(Chem.MolToSmiles(h))].append(s)
print(dict(grp))
assert sorted(len(v) for v in grp.values()) == [3, 3]
assert canon('CCCCC') in grp and canon('CC(C)CC') in grp and canon('CC(C)(C)C') not in grp
c3 = enumerate_isomers({'C': 3, 'Cl': 1}, 5)
c3st = [x for s in c3 for x in stereo_isomers(s)]
print('C3H5Cl:', len(c3), len(c3st), c3st); assert len(c3st) == 5
c4 = enumerate_isomers({'C': 4, 'Cl': 1}, 9)
c4st = [x for s in c4 for x in stereo_isomers(s)]
print('C4H9Cl:', len(c4), len(c4st)); assert (len(c4), len(c4st)) == (4, 5)
print('ORG-STR-007 OK')
