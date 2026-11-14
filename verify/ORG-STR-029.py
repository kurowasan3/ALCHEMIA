# 検証スクリプト ORG-STR-029  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/ORG-STR-029/item.json'
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
def all_equations():
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
def mirror(s): return can(s.replace('@@', '!').replace('@', '@@').replace('!', '@'))
acids = {}
for s in ['CCCC', 'CC(C)C']:
    m = Chem.MolFromSmiles(s); n = m.GetNumAtoms()
    for a in range(n):
        for b in range(a, n):
            rw = Chem.RWMol(m)
            for i in (a, b):
                c = rw.AddAtom(Chem.Atom(6)); rw.AddBond(i, c, Chem.BondType.SINGLE)
                o1 = rw.AddAtom(Chem.Atom(8)); rw.AddBond(c, o1, Chem.BondType.DOUBLE)
                o2 = rw.AddAtom(Chem.Atom(8)); rw.AddBond(c, o2, Chem.BondType.SINGLE)
            mm = rw.GetMol()
            try: Chem.SanitizeMol(mm)
            except Exception: continue
            acids[Chem.MolToSmiles(mm)] = mm
assert len(acids) == 9
achiral = []
for k, m in acids.items():
    for st in set(Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m)):
        if mirror(st) == can(st): achiral.append(st)
assert len(achiral) == 7, achiral
def monoesters(st):
    mm = Chem.MolFromSmiles(st); out = set()
    for mt in mm.GetSubstructMatches(Chem.MolFromSmarts('C(=O)[OX2H1]')):
        rw = Chem.RWMol(mm); c = rw.AddAtom(Chem.Atom(6)); rw.AddBond(mt[2], c, Chem.BondType.SINGLE)
        x = rw.GetMol(); Chem.SanitizeMol(x); out.add(Chem.MolToSmiles(x))
    return out
chiral_mono = []
for st in achiral:
    ok = False
    for mo in monoesters(st):
        centers = Chem.FindMolChiralCenters(Chem.MolFromSmiles(mo), includeUnassigned=True, useLegacyImplementation=False)
        if any(c[1] == '?' for c in centers) or mirror(mo) != can(mo): ok = True
    if ok: chiral_mono.append(st)
assert len(chiral_mono) == 5, chiral_mono
assert can('O=C(O)CCCCC(=O)O') not in chiral_mono and can('CC(C)(CC(=O)O)C(=O)O') not in chiral_mono
meso = can('C[C@H](C(=O)O)[C@@H](C)C(=O)O'); assert meso in achiral and meso in chiral_mono
ms = sorted(monoesters(meso)); assert len(ms) == 2 and mirror(ms[0]) == ms[1]
E = can('C[C@H](C(=O)OC)[C@@H](C)C(=O)OC'); assert mirror(E) == E
# 酸無水物の環の員数
for smi, size in [('CC1C(C)C(=O)OC1=O', 5), ('CC1CC(=O)OC(=O)C1', 6)]:
    assert [len(r) for r in Chem.GetSymmSSSR(Chem.MolFromSmiles(smi))] == [size]
assert check_equation('CH3CH2CH2CH(COOH)2 -> CH3CH2CH2CH2COOH + CO2') is None
assert check_equation('C7H12O4 + H2O -> C6H10O4 + CH4O') is None
print('対称性の破れ・メソ形 OK')

print('ALL OK')
