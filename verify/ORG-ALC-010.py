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

# ORG-ALC-010 検証：C5H12O アルコールの分子内脱水で生じるアルケン（立体異性体を区別）
alcohols = [s for s in enumerate_isomers({'C': 5, 'O': 1}, 12) if has(s, '[CX4][OX2H]')]
assert len(alcohols) == 8
def dehydrate(smi):
    m = Chem.MolFromSmiles(smi)
    prods = set()
    for (c, o) in m.GetSubstructMatches(Chem.MolFromSmarts('[CX4:1]-[OX2H:2]')):
        for nb in m.GetAtomWithIdx(c).GetNeighbors():
            if nb.GetSymbol() == 'C' and nb.GetTotalNumHs() > 0:
                rw = Chem.RWMol(m)
                rw.GetBondBetweenAtoms(c, nb.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
                rw.RemoveAtom(o)
                for a in rw.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
                mm = rw.GetMol(); Chem.SanitizeMol(mm)
                for x in stereo_isomers(Chem.MolToSmiles(mm)):
                    prods.add(x)
    return sorted(prods)
res = {s: dehydrate(s) for s in alcohols}
for k, v in res.items(): print(k, '->', v)
assert res[canon('CC(O)CCC')] == sorted([canon('C=CCCC'), 'C/C=C/CC', 'C/C=C\\CC'])
assert res[canon('CCC(O)CC')] == sorted(['C/C=C/CC', 'C/C=C\\CC'])
assert sorted(res[canon('CCC(C)(C)O')]) == sorted([canon('C=C(C)CC'), canon('CC=C(C)C')])
none = [k for k, v in res.items() if not v]; assert none == [canon('CC(C)(C)CO')]
one = [k for k, v in res.items() if len(v) == 1]; print('one alkene:', one); assert len(one) == 3
# ザイツェフ則：二重結合の炭素についた H の数が少ない方が主
def hs_on_db(smi):
    m = Chem.MolFromSmiles(smi); b = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE][0]
    return b.GetBeginAtom().GetTotalNumHs() + b.GetEndAtom().GetTotalNumHs()
major = min(res[canon('CCC(C)(C)O')], key=hs_on_db); assert major == canon('CC=C(C)C')
print('ORG-ALC-010 OK')
