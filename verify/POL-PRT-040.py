# POL-PRT-040 検証：グリシン2個とリシン1個からなる鎖状トリペプチド（側鎖のε-アミノ基による結合も含む）の数え上げ
import itertools
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# 残基: (SMILES, {'N': [アミノ基 N のidx...], 'C': [(カルボニルC, OH の O)...]})
RES = {
 'Gly': ('NCC(=O)O', {'N': {'a': 0}, 'C': {'a': (2, 4)}}),
 'Lys': ('NC(CCCCN)C(=O)O', {'N': {'a': 0, 'e': 6}, 'C': {'a': (7, 9)}}),
 'Ala': ('NC(C)C(=O)O', {'N': {'a': 0}, 'C': {'a': (3, 5)}}),
}

def enumerate_peptides(names):
    mols = [Chem.MolFromSmiles(RES[n][0]) for n in names]
    combo = Chem.MolFromSmiles('.'.join(RES[n][0] for n in names))
    off = []; o = 0
    for m in mols: off.append(o); o += m.GetNumAtoms()
    cgroups = [(i, k, off[i] + c, off[i] + oh) for i, n in enumerate(names) for k, (c, oh) in RES[n][1]['C'].items()]
    ngroups = [(i, k, off[i] + a) for i, n in enumerate(names) for k, a in RES[n][1]['N'].items()]
    edges = [(c, n) for c in cgroups for n in ngroups if c[0] != n[0]]
    found = {}
    for es in itertools.combinations(edges, len(names) - 1):
        cs = [e[0][:2] for e in es]; ns = [e[1][:2] for e in es]
        if len(set(cs)) < len(cs) or len(set(ns)) < len(ns): continue
        # 連結（木）か
        parent = list(range(len(names)))
        def f(x):
            while parent[x] != x: x = parent[x]
            return x
        ok = True
        for e in es:
            a, b = f(e[0][0]), f(e[1][0])
            if a == b: ok = False; break
            parent[a] = b
        if not ok: continue
        rw = Chem.RWMol(combo)
        rm = []
        for (ci, ck, c, oh), (ni, nk, n) in es:
            rw.AddBond(c, n, Chem.BondType.SINGLE); rm.append(oh)
        for a in sorted(rm, reverse=True): rw.RemoveAtom(a)
        m = rw.GetMol(); Chem.SanitizeMol(m)
        smi = Chem.MolToSmiles(m)
        found.setdefault(smi, []).append([(names[e[0][0]] + e[0][1], names[e[1][0]] + e[1][1]) for e in es])
    return found

if __name__ == '__main__':
    gk = enumerate_peptides(['Gly', 'Gly', 'Lys'])
    print('Gly2Lys 構造異性体数', len(gk))
    free_alpha = 0
    for s, how in gk.items():
        m = Chem.MolFromSmiles(s)
        print(' ', s, how[0], rdMolDescriptors.CalcMolFormula(m))
    # Lys の α-アミノ基が遊離しているもの: Lys a の N が結合に使われていない
    fa = [s for s, how in gk.items() if not any(b == 'Lysa' for _, b in how[0])]
    eps = [s for s, how in gk.items() if any(b == 'Lyse' for _, b in how[0])]
    print('Lys α-NH2 遊離', len(fa), ' ε結合を含む', len(eps))
    assert len(gk) == 6 and len(fa) == 3 and len(eps) == 3
    ga = enumerate_peptides(['Gly', 'Gly', 'Ala'])
    print('Gly2Ala', len(ga)); assert len(ga) == 3
    # 遊離の官能基: 全て同じ分子式 C10H20N4O4 ，遊離 NH2 2個・COOH 1個
    for s in gk:
        m = Chem.MolFromSmiles(s)
        nh2 = len(m.GetSubstructMatches(Chem.MolFromSmarts('[NX3;H2][CX4]')))
        cooh = len(m.GetSubstructMatches(Chem.MolFromSmarts('C(=O)[OH]')))
        assert (nh2, cooh) == (2, 1), (s, nh2, cooh)
        assert rdMolDescriptors.CalcMolFormula(m) == 'C10H20N4O4'
    # 不斉炭素は各1個（Lys の α 炭素）→ 立体異性体を含めると 12
    from rdkit.Chem import FindMolChiralCenters
    for s in gk:
        assert len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False)) == 1
    print('分子量', round(Descriptors.MolWt(Chem.MolFromSmiles(list(gk)[0])), 1))
    print('ALL OK')
