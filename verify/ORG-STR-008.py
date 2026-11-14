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

# ORG-STR-008 検証：C4H8O2 の異性体から条件に合うものを選ぶ（一意性の確認）
allc = enumerate_isomers({'C': 4, 'O': 2}, 8)
bad = ['OO', '[CX3]=[CX3][OX2H]', '[CX4]([OX2H])[OX2H]', '[CX4]([OX2H])[OX2]', 'C=C=C', 'C=C=O', '[#6]=[#6]-[OX2]-[OX2]']
S = [s for s in allc if not any(has(s, b) for b in bad)]
print(len(allc), '->', len(S), 'after excluding unstable')
ring = Chem.MolFromSmarts('[R]')
chain = [s for s in S if not Chem.MolFromSmiles(s).HasSubstructMatch(ring)]
COOH = '[CX3](=O)[OX2H]'; ESTER = '[#6][CX3](=O)[OX2][#6]'; FORMATE = '[CX3H1](=O)[OX2][#6]'
CHO = '[CX3H1](=O)[#6]'; OH = '[CX4][OX2H]'
def branched(s):
    m = Chem.MolFromSmiles(s); return any(sum(1 for n in a.GetNeighbors() if n.GetSymbol() == 'C') >= 3 for a in m.GetAtoms() if a.GetSymbol() == 'C')
def hydrolyze(s):
    m = Chem.MolFromSmiles(s); match = m.GetSubstructMatch(Chem.MolFromSmarts('[CX3:1](=O)[OX2:2][#6:3]'))
    c, _, o, r = match
    rw = Chem.RWMol(m); rw.RemoveBond(c, o); n = rw.AddAtom(Chem.Atom(8)); rw.AddBond(c, n, Chem.BondType.SINGLE)
    for a in rw.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    mm = rw.GetMol(); Chem.SanitizeMol(mm)
    fr = [Chem.MolToSmiles(f) for f in Chem.GetMolFrags(mm, asMols=True)]
    acid = [f for f in fr if has(f, '[CX3](=O)[OX2H]')][0]; alc = [f for f in fr if f != acid][0]
    return acid, alc
silver = lambda s: has(s, '[CX3H1]=O')
iodo_alc = lambda s: has(s, '[CH3][CX4H1,CX4H2][OH]') and not has(s, '[CH3][CH2][CH2]')  # CH3CH(OH)- のみ（エタノール含む）
iodo = lambda s: has(s, '[CH3][CX4H1]([OH])') or has(s, '[CH3][CX3](=O)[#6,#1]')
secondary = lambda s: has(s, '[CX4H1]([#6])([#6])[OX2H]')
r1 = [s for s in S if has(s, COOH) and branched(s)]
esters = [s for s in S if has(s, ESTER) or has(s, FORMATE)]
# ヨードホルムを示すアルコール：エタノール，2-プロパノール（C4H8O2 のエステルから出うるもの）
r2 = [s for s in esters if (lambda a, b: not silver(a) and canon(b) in {canon('CCO'), canon('CC(C)O')})(*hydrolyze(s))]
r3 = [s for s in esters if (lambda a, b: silver(a) and secondary(b))(*hydrolyze(s))]
hc = [s for s in chain if len(Chem.MolFromSmiles(s).GetSubstructMatches(Chem.MolFromSmarts(OH))) == 1 and len(Chem.MolFromSmiles(s).GetSubstructMatches(Chem.MolFromSmarts(CHO))) + has(s, '[CH2]=O') == 1 and not has(s, 'C=C')]
print('hydroxy aldehydes:', hc)
r4 = [s for s in hc if n_chiral(s) > 0 and has(s, '[CH3][CX4H1][OH]')]
r5 = [s for s in hc if branched(s) and n_chiral(s) == 0]
diols = [s for s in chain if len(Chem.MolFromSmiles(s).GetSubstructMatches(Chem.MolFromSmarts(OH))) == 2]
r6 = [s for s in diols if len(stereo_isomers(s)) == 2 and any('/' in x for x in stereo_isomers(s))]
print('diols:', diols)
res = [r1, r2, r3, r4, r5, r6]; print(res)
exp = ['CC(C)C(=O)O', 'CCOC(C)=O', 'CC(C)OC=O', 'CC(O)CC=O', 'CC(C)(O)C=O', 'OCC=CCO']
for r, e in zip(res, exp):
    assert r == [canon(e)], (r, e)
trans = [x for x in stereo_isomers('OCC=CCO') if '/' in x]; print('stereo:', trans)
print('ORG-STR-008 OK')
