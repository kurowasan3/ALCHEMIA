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

# ORG-ALC-011 検証：C5H12O アルコールの不斉炭素と，酸化・脱水＋水素付加での変化
alc = sorted(s for s in enumerate_isomers({'C': 5, 'O': 1}, 12) if has(s, '[CX4][OX2H]')); assert len(alc) == 8
chiral = [s for s in alc if n_chiral(s) > 0]; print('chiral:', chiral)
assert sorted(chiral) == sorted(canon(x) for x in ['CC(O)CCC', 'OCC(C)CC', 'CC(O)C(C)C'])
def oxidize(s):
    m = Chem.MolFromSmiles(s)
    c, o = m.GetSubstructMatch(Chem.MolFromSmarts('[CX4:1][OX2H]'))
    if m.GetAtomWithIdx(c).GetTotalNumHs() == 0: return None   # 第三級は酸化されない
    rw = Chem.RWMol(m); rw.GetBondBetweenAtoms(c, o).SetBondType(Chem.BondType.DOUBLE)
    for a in rw.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    mm = rw.GetMol(); Chem.SanitizeMol(mm); return Chem.MolToSmiles(mm)
for s in chiral:
    p = oxidize(s); print(s, '->', p, 'chiral' if n_chiral(p) else 'achiral')
kept = [s for s in chiral if n_chiral(oxidize(s)) > 0]; assert kept == [canon('OCC(C)CC')]
notox = [s for s in alc if oxidize(s) is None]; assert notox == [canon('CCC(C)(C)O')]
# P：不斉炭素あり・ヨードホルム陽性（CH3CH(OH)-）・脱水→水素付加で枝分かれアルカン
iodo = Chem.MolFromSmarts('[CH3][CH1]([OH])[#6]')
def skeleton(s):
    m = Chem.RWMol(Chem.MolFromSmiles(s)); o = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    m.RemoveAtom(o); mm = m.GetMol(); Chem.SanitizeMol(mm); return canon(Chem.MolToSmiles(mm))
P = [s for s in chiral if Chem.MolFromSmiles(s).HasSubstructMatch(iodo) and skeleton(s) != canon('CCCCC')]
print('P', P); assert P == [canon('CC(O)C(C)C')]
print('ORG-ALC-011 OK')
