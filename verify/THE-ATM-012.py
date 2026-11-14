# 検証スクリプト THE-ATM-012  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-ATM-012/item.json'
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
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b:
                eqs.add(b)
    return sorted(eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
import itertools
from sympy import Rational as R, symbols, simplify
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
check_structure()
# D を含む式：D→H に置きかえて原子・電荷を確認し，D の数も左右で数える
def countD(side):
    n = 0
    for term in side.split(' + '):
        term = term.strip(); m = re.match(r'(\d*)(.*)', term); k = int(m.group(1) or 1); f = m.group(2)
        n += k * sum(int(x or 1) for x in re.findall(r'D(\d*)', f))
    return n
eqs = all_equations(); print(eqs)
assert len(eqs) == 3
for e in eqs:
    assert check_equation(e.replace('D', 'H')) is None, e
    l, r_ = re.split(r'<=>|->', e); assert countD(l) == countD(r_), e
# (1)
f1 = R(360, 100) / (R(360,100) + R(20,100)); assert f1 == R(18, 19)
M1 = R(355,10) + 1*(1-f1) + 2*f1; print('M1', float(M1)); assert sig(float(M1), 3) == 37.4
# (2)
f = symbols('f', positive=True)
K = simplify((2*f*(1-f))**2 / ((1-f)**2 * f**2)); assert K == 4
nHDO = R(18,10) * 2*f1*(1-f1); print('HDO', float(nHDO)); assert sig(float(nHDO), 2) == 0.18
# (3)
f3 = R(2) / (R(1,2) + 2); assert f3 == R(4,5)
M3 = 60 + f3; assert sig(float(M3), 3) == 60.8
# (4)
Mx = 3*12 + 8 + 2*16; assert Mx == 76
vals = {}
for n in range(0, 3):
    fn = R(8,10) / (R(8,10) + R(1,10)*n); vals[n] = Mx + n*fn
print({k: float(v) for k, v in vals.items()})
assert [n for n, v in vals.items() if sig(float(v), 3) == 77.6] == [2]
# C3H8O2 の構造列挙（重原子 C,C,C,O,O の木）
labels = ['C', 'C', 'C', 'O', 'O']
smis = set()
for prufer in itertools.product(range(5), repeat=3):
    deg = [1]*5
    for x in prufer: deg[x] += 1
    seq = list(prufer); edges = []; d = deg[:]
    for x in seq:
        leaf = min(i for i in range(5) if d[i] == 1)
        edges.append((leaf, x)); d[leaf] -= 1; d[x] -= 1
    u, v = [i for i in range(5) if d[i] == 1]; edges.append((u, v))
    rw = Chem.RWMol()
    for l in labels: rw.AddAtom(Chem.Atom(l))
    for a, b in edges: rw.AddBond(a, b, Chem.BondType.SINGLE)
    m = rw.GetMol()
    try:
        Chem.SanitizeMol(m)
    except Exception:
        continue
    smis.add(Chem.MolToSmiles(m))
mols = {s: Chem.MolFromSmiles(s) for s in smis}
for s, m in mols.items(): assert rdMolDescriptors.CalcMolFormula(m) == 'C3H8O2'
print('全異性体（O–O 含む）', len(smis), sorted(smis))
def ok(m):
    for b in m.GetBonds():
        if b.GetBeginAtom().GetSymbol() == 'O' and b.GetEndAtom().GetSymbol() == 'O': return False
    for a in m.GetAtoms():
        if a.GetSymbol() == 'C' and sum(1 for n in a.GetNeighbors() if n.GetSymbol() == 'O') >= 2: return False
    return True
cand = sorted(s for s, m in mols.items() if ok(m)); print('条件を満たす', cand)
assert len(cand) == 3
patt = Chem.MolFromSmarts('[OX2H]')
def nOH(m): return len(m.GetSubstructMatches(patt))
def chiral(m): return len(Chem.FindMolChiralCenters(m, includeUnassigned=True)) > 0
sel = [s for s in cand if nOH(mols[s]) == 2 and chiral(mols[s])]; print('X', sel)
assert sel == [Chem.MolToSmiles(Chem.MolFromSmiles('CC(O)CO'))]
assert sorted(cand) == sorted(Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in ['CC(O)CO', 'OCCCO', 'COCCO'])
print('全検証OK:', ITEM)
