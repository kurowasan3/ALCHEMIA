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

# ORG-HYD-011 検証：C6H12 鎖式アルケンの KMnO4 酸化生成物と水素付加生成物から P, Q, R を一意に決める
from sympy import symbols, solve
n = symbols('n'); assert solve(n + 2 * n / 4 - 9, n) == [6]       # CnH2n + 3n/2 O2
ring = Chem.MolFromSmarts('[R]')
alk = [s for s in enumerate_isomers({'C': 6}, 12) if not Chem.MolFromSmiles(s).HasSubstructMatch(ring)]
print(len(alk), 'C6H12 alkenes (constitutional)')
def cleave(smi):
    m = Chem.MolFromSmiles(smi)
    db = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE][0]
    prods = []
    rw = Chem.RWMol(m); a, b = db.GetBeginAtomIdx(), db.GetEndAtomIdx()
    rw.RemoveBond(a, b)
    frags = []
    for idx in (a, b):
        r = Chem.RWMol(rw)
        o = r.AddAtom(Chem.Atom(8)); r.AddBond(idx, o, Chem.BondType.DOUBLE)
        mm = r.GetMol(); Chem.SanitizeMol(mm)
        frag = [f for f in Chem.GetMolFrags(mm, asMols=True, sanitizeFrags=True) if f.GetNumAtoms() and any(at.GetSymbol()=='O' for at in f.GetAtoms())]
        # choose fragment containing idx
        for f in Chem.GetMolFrags(mm, asMols=False):
            if idx in f:
                sub = Chem.RWMol(mm)
                for k in sorted(set(range(mm.GetNumAtoms())) - set(f), reverse=True): sub.RemoveAtom(k)
                s = Chem.MolToSmiles(sub.GetMol())
        hs = m.GetAtomWithIdx(idx).GetTotalNumHs()
        if hs == 2: frags.append('CO2')
        elif hs == 1:
            s2 = Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('C=O', 'C(=O)O', 1))) if False else None
            mm2 = Chem.RWMol(Chem.MolFromSmiles(s))
            # aldehyde -> acid
            patt = Chem.MolFromSmarts('[CX3H1]=O'); match = mm2.GetSubstructMatch(patt)
            oo = mm2.AddAtom(Chem.Atom(8)); mm2.AddBond(match[0], oo, Chem.BondType.SINGLE)
            frags.append(('acid', Chem.MolToSmiles(mm2.GetMol())))
        else:
            frags.append(('ketone', s))
    return frags
def hydro(smi):
    h = Chem.RWMol(Chem.MolFromSmiles(smi))
    for b in h.GetBonds(): b.SetBondType(Chem.BondType.SINGLE)
    return canon(Chem.MolToSmiles(h))
def n_monochloro(alkane):
    m = Chem.MolFromSmiles(alkane)
    ranks = Chem.CanonicalRankAtoms(m, breakTies=False)
    return len({ranks[a.GetIdx()] for a in m.GetAtoms() if a.GetTotalNumHs() > 0})
info = {s: (cleave(s), hydro(s)) for s in alk}
def kinds(fr): return sorted(x if x == 'CO2' else x[0] for x in fr)
P = [s for s, (fr, h) in info.items() if kinds(fr) == ['ketone', 'ketone'] and len({x[1] for x in fr}) == 1]
Q = [s for s, (fr, h) in info.items() if kinds(fr) == ['CO2', 'ketone'] and h in {info[p][1] for p in P}]
Rc = [s for s, (fr, h) in info.items() if kinds(fr) == ['CO2', 'acid']]
R = [s for s in Rc if info[s][1] not in {info[p][1] for p in P} and n_monochloro(info[s][1]) == 3 and info[s][1] != canon('CCCCCC')]
print('P', P, 'Q', Q, 'R candidates', Rc, '-> R', R)
assert P == [canon('CC(C)=C(C)C')] and Q == [canon('C=C(C)C(C)C')] and R == [canon('C=CC(C)(C)C')]
print('S', info[P[0]][1], 'T', info[R[0]][1], 'monochloro(T)=', n_monochloro(info[R[0]][1]))
print('Q products', info[Q[0]][0])
# ヨードホルム：CH3-CO- をもつか
ket = [x[1] for x in info[Q[0]][0] if x != 'CO2'][0]
assert Chem.MolFromSmiles(ket).HasSubstructMatch(Chem.MolFromSmarts('[CH3]C(=O)[#6]'))
print('ORG-HYD-011 OK')
