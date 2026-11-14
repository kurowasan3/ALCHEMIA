ITEM = '/home/claude/chem2/outh/items/ORG-STR-021/item.json'
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
from rdkit import Chem, RDLogger
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
RDLogger.DisableLog('rdApp.*')
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
C4H8 = ['C=CCC', 'CC=CC', 'C=C(C)C', 'C1CCC1', 'CC1CC1']
def subst(s):
    m = Chem.AddHs(Chem.MolFromSmiles(s)); out = set()
    for a in m.GetAtoms():
        if a.GetSymbol() == 'H':
            rw = Chem.RWMol(m); rw.GetAtomWithIdx(a.GetIdx()).SetAtomicNum(17)
            out.add(Chem.MolToSmiles(Chem.RemoveHs(rw)))
    return out
iso = set().union(*[subst(s) for s in C4H8])
ring = {s for s in iso if Chem.MolFromSmiles(s).GetRingInfo().NumRings()}
print(len(iso), len(iso - ring), len(ring)); assert (len(iso - ring), len(ring)) == (8, 4)
def stereo(s): return sorted({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(s), options=StereoEnumerationOptions(unique=True))})
mirror = lambda s: can(s.replace('@@', '!').replace('@', '@@').replace('!', '@'))
tot = sum(len(stereo(s)) for s in iso); print('total', tot); assert tot == 19
chiral = sorted(s for s in iso if any(mirror(x) != x for x in stereo(s)))
print('chiral', chiral); assert chiral == sorted([can('C=CC(C)Cl'), can('CC1CC1Cl')])
assert len(stereo('CC1CC1Cl')) == 4
dm = stereo('CC1CC1C'); assert len(dm) == 3 and sum(mirror(x) == x for x in dm) == 1  # メソ 1 つ
def cistrans(s): return any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetStereo() != Chem.BondStereo.STEREONONE for b in Chem.MolFromSmiles(stereo(s)[0]).GetBonds())
def hyd(s):
    m = Chem.RWMol(Chem.MolFromSmiles(s))
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE: b.SetBondType(Chem.BondType.SINGLE)
    for a in m.GetAtoms(): a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    mm = m.GetMol(); Chem.SanitizeMol(mm); return Chem.MolToSmiles(mm)
def has_chiral_C(s):
    return bool(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
ct = sorted(s for s in iso if cistrans(s)); print('cis-trans', ct); assert len(ct) == 3
X = [s for s in ct if s not in chiral and has_chiral_C(hyd(s))]; print('X', X); assert X == [can('CC=C(C)Cl')]
prods = {hyd(s) for s in iso - ring}; print(prods)
assert can('CC(C)(C)Cl') not in prods
assert sorted(s for s in iso - ring if hyd(s) == can('CC(C)CCl')) == sorted([can('CC(C)=CCl'), can('C=C(C)CCl')])

print('全検証OK:', ITEM)
