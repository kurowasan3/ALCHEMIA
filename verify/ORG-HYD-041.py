# ORG-HYD-041 検証：C6H12・C7H14 アルケンへの Br2・H2 付加で生じる不斉炭素の数
import sys; sys.path.insert(0, "/home/claude/kurowasan3/alchemia/tools")
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')
from itertools import combinations
import sys

def grow(s):
    out=set()
    for x in s:
        m=Chem.MolFromSmiles(x)
        for a in m.GetAtoms():
            rw=Chem.RWMol(m); i=rw.AddAtom(Chem.Atom(6)); rw.AddBond(a.GetIdx(),i,Chem.BondType.SINGLE)
            mm=rw.GetMol()
            try: Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
            except: pass
    return out

def nch(s): return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s),includeUnassigned=True,useLegacyImplementation=True))
def add(s, X):
    m=Chem.MolFromSmiles(s)
    rxn=AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1](%s)[C:2](%s)'%(X,X))
    p=rxn.RunReactants((m,))[0][0]; Chem.SanitizeMol(p); return Chem.MolToSmiles(p)
def alkenes(n):
    lv = {'C'}
    for k in range(n - 1): lv = grow(lv)
    out = set()
    for s in lv:
        m = Chem.MolFromSmiles(s)
        for b in m.GetBonds():
            rw = Chem.RWMol(m); rw.GetBondBetweenAtoms(b.GetBeginAtomIdx(), b.GetEndAtomIdx()).SetBondType(Chem.BondType.DOUBLE)
            mm = rw.GetMol()
            try: Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
            except Exception: pass
    return out, lv
a6, alk6 = alkenes(6)
rows = []
for s in sorted(a6):
    if nch(s) > 0: continue
    rows.append((s, nch(add(s, 'Br')), nch(Chem.MolToSmiles(Chem.MolFromSmiles(add(s, '[H]'))))))
for r in rows: print(r)
two = sorted(r[0] for r in rows if r[1] == 2); zero = sorted(r[0] for r in rows if r[1] == 0)
print('(1)', two, '(2)', zero)
assert two == sorted(['CC=C(C)CC', 'CC=CC(C)C', 'CC=CCCC', 'CCC=CCC']) and zero == sorted(['C=C(CC)CC', 'CC(C)=C(C)C'])
assert all(r[2] == 0 for r in rows)
print('C6H14 のアルカンの不斉炭素', [(s, nch(s)) for s in sorted(alk6)]); assert all(nch(s) == 0 for s in alk6)
a7, alk7 = alkenes(7)
print('C7H16 で不斉炭素をもつもの', [s for s in sorted(alk7) if nch(s) > 0])
W = []
for s in sorted(a7):
    if nch(s) > 0: continue
    m = Chem.MolFromSmiles(s)
    noH = all(a.GetTotalNumHs() == 0 for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE for a in (b.GetBeginAtom(), b.GetEndAtom()))
    h = Chem.MolToSmiles(Chem.MolFromSmiles(add(s, '[H]')))
    if nch(h) > 0 and noH: W.append((s, nch(add(s, 'Br')), h))
print('(4)', W); assert W == [('CCC(C)=C(C)C', 1, 'CCC(C)C(C)C')]
from chemlib import check_equation
for e in ['C6H12 + Br2 -> C6H12Br2','C7H14 + H2 -> C7H16']: assert check_equation(e) is None
print('OK')
