# ORG-HYD-040 検証：C5H10 アルケンへの HCl 付加（マルコフニコフ則）と不斉炭素
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from rdkit import Chem
for n in range(2, 8):
    print(n, round((14*n+36.5)/(14*n), 4))
assert [n for n in range(2, 9) if round((14*n+36.5)/(14*n), 2) == 1.52] == [5]
alk = {'1-ペンテン': 'C=CCCC', '2-ペンテン': 'CC=CCC', '2-メチル-1-ブテン': 'C=C(C)CC',
       '3-メチル-1-ブテン': 'C=CC(C)C', '2-メチル-2-ブテン': 'CC=C(C)C'}
def markov(smi):
    m = Chem.MolFromSmiles(smi); b = [x for x in m.GetBonds() if x.GetBondTypeAsDouble() == 2][0]
    a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
    h1, h2 = a1.GetTotalNumHs(), a2.GetTotalNumHs()
    targets = [a1.GetIdx()] if h1 < h2 else [a2.GetIdx()] if h2 < h1 else [a1.GetIdx(), a2.GetIdx()]  # Cl は H の少ない側
    prods = set()
    for t in targets:
        em = Chem.RWMol(m); em.GetBondBetweenAtoms(a1.GetIdx(), a2.GetIdx()).SetBondType(Chem.BondType.SINGLE)
        cl = em.AddAtom(Chem.Atom(17)); em.AddBond(t, cl, Chem.BondType.SINGLE)
        p = em.GetMol(); Chem.SanitizeMol(p); prods.add(Chem.MolToSmiles(p))
    return prods
res = {}
for k, s in alk.items():
    ps = markov(s)
    ch = {p: len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(p), includeUnassigned=True, useLegacyImplementation=False)) for p in ps}
    res[k] = ch; print(k, ch)
ok = [k for k, ch in res.items() if len(ch) == 1 and list(ch.values())[0] == 0]
print('single achiral product:', ok)
assert ok == ['2-メチル-1-ブテン', '2-メチル-2-ブテン']
# C=C 炭素上の H の数
for k in ok:
    m = Chem.MolFromSmiles(alk[k]); b = [x for x in m.GetBonds() if x.GetBondTypeAsDouble() == 2][0]
    print(k, 'H on C=C', b.GetBeginAtom().GetTotalNumHs() + b.GetEndAtom().GetTotalNumHs())
for e in ['CH3CH=C(CH3)2 + HCl -> CH3CH2CCl(CH3)2']:
    assert check_equation(e) is None
print('OK')
