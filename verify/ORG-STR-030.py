# 検証スクリプト ORG-STR-030  python3 -I この.py
ITEM = '/home/claude/chem3/outh/items/ORG-STR-030/item.json'
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
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
mirror = lambda s: can(s.replace('@@', '!').replace('@', '@@').replace('!', '@'))
W = {'C': 12, 'H': 1.0, 'O': 16}
Mw = lambda s: sum(W[a.GetSymbol()] for a in Chem.AddHs(Chem.MolFromSmiles(s)).GetAtoms())
A = 'O1C(C)C(C)OC(=O)CCCCC(=O)OC(C)CC1=O'
assert CalcMolFormula(Chem.MolFromSmiles(A)) == 'C14H22O6' and Mw(A) == 286
assert [len(r) for r in Chem.GetSymmSSSR(Chem.MolFromSmiles(A))] == [14]
assert Mw('CC(O)C(C)O') == 90 and Mw('CC(O)CC(=O)O') == 104 and Mw('OC(=O)CCCCC(=O)O') == 146
# D の候補と環の員数
acids = set()
for s in ['CCCC', 'CC(C)C']:
    m = Chem.MolFromSmiles(s); n = m.GetNumAtoms()
    for a in range(n):
        for b in range(a, n):
            rw = Chem.RWMol(m)
            for i in (a, b):
                c = rw.AddAtom(Chem.Atom(6)); rw.AddBond(i, c, Chem.BondType.SINGLE)
                o1 = rw.AddAtom(Chem.Atom(8)); rw.AddBond(c, o1, Chem.BondType.DOUBLE)
                o2 = rw.AddAtom(Chem.Atom(8)); rw.AddBond(c, o2, Chem.BondType.SINGLE)
            x = rw.GetMol()
            try: Chem.SanitizeMol(x)
            except Exception: continue
            acids.add(Chem.MolToSmiles(x))
assert len(acids) == 9
ring = {}
for s in acids:
    m = Chem.MolFromSmiles(s)
    cs = [mt[0] for mt in m.GetSubstructMatches(Chem.MolFromSmarts('[CX3](=O)[OX2H1]'))]
    k = len(Chem.GetShortestPath(m, cs[0], cs[1])) - 2
    ring[s] = k + 10
assert [s for s, r in ring.items() if r == 14] == [can('OC(=O)CCCCC(=O)O')]
# C の候補
hyd = ['CCC(O)C(=O)O', 'CC(O)CC(=O)O', 'OCCCC(=O)O', 'CC(C)(O)C(=O)O', 'OCC(C)C(=O)O']
nch = lambda s: len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
iod = Chem.MolFromSmarts('[CH3][CX4H1][OX2H1]')
assert [s for s in hyd if nch(s) and Chem.MolFromSmiles(s).HasSubstructMatch(iod)] == ['CC(O)CC(=O)O']
# 立体異性体
def hydrolyze(smi):
    m = Chem.MolFromSmiles(smi); rw = Chem.RWMol(m)
    for a, _, b, _ in m.GetSubstructMatches(Chem.MolFromSmarts('[CX3](=O)-[OX2]-[#6]')):
        rw.RemoveBond(a, b); o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(a, o, Chem.BondType.SINGLE)
    x = rw.GetMol(); Chem.SanitizeMol(x)
    return [can(f) for f in Chem.MolToSmiles(x).split('.')]
iso = sorted(set(Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(A))))
assert len(iso) == 8
def Bpart(s): return [f for f in hydrolyze(s) if can(Chem.MolToSmiles(Chem.MolFromSmiles(f))) and f.count('O') == 2 and 'C(=O)' not in f][0]
def Cpart(s): return [f for f in hydrolyze(s) if f.count('O') == 3][0]
meso = [s for s in iso if mirror(Bpart(s)) == Bpart(s)]
assert len(meso) == 4 and all(mirror(s) in meso and mirror(s) != s for s in meso)
Cs = sorted({Cpart(s) for s in meso}); assert len(Cs) == 2 and mirror(Cs[0]) == Cs[1]
fixed = [s for s in meso if Cpart(s) == Cs[0]]
assert len(fixed) == 2 and mirror(fixed[0]) != fixed[1]
assert check_equation('C14H22O6 + 3H2O -> C4H10O2 + C4H8O3 + C6H10O4') is None
print('環の員数・成分・立体異性体 OK')

print('ALL OK')
