# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out2/items/ORG-STR-013/item.json'
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and not re.search(r'_n|\bn[A-Z]|R', b):
            eqs.add(b)
def iso_count(side):
    # 同位体 ^{18}O などの数（係数込み）を数える
    tot = 0
    for term in side.split(' + '):
        term = term.strip(); mm = re.match(r'^(\d+)', term); k = int(mm.group(1)) if mm else 1
        tot += k * len(re.findall(r'\^\{\d+\}', term))
    return tot
for e in sorted(eqs):
    plain = re.sub(r'\^\{\d+\}(?=[A-Z])', '', e)
    r = check_equation(plain)
    assert r is None, (e, r)
    if plain != e:
        L, Rr = re.split(r'->|<=>', e)
        assert iso_count(L) == iso_count(Rr), ('同位体の数', e)
    print('式OK:', e)

# 小問数と答えの数の一致
nq = it['problem'].count('\\item')
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
print('小問', nq, '答え', na); assert nq == na or (nq == 0 and na == 0)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

# ---- 立体化学の補助（RDKit） ----
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
RDLogger = __import__('rdkit.RDLogger', fromlist=['RDLogger']); RDLogger.DisableLog('rdApp.*')

def can(m):
    return Chem.MolToSmiles(m, isomericSmiles=True)

def mirror(smi):
    m = Chem.MolFromSmiles(smi)
    for a in m.GetAtoms():
        t = a.GetChiralTag()
        if t == Chem.ChiralType.CHI_TETRAHEDRAL_CW: a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW: a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    return can(m)

def stereoisomers(smi):
    m = Chem.MolFromSmiles(smi)
    opts = StereoEnumerationOptions(unique=True, onlyUnassigned=True)
    return sorted({can(x) for x in EnumerateStereoisomers(m, options=opts)})

def is_chiral(smi):
    s = can(Chem.MolFromSmiles(smi))
    return mirror(s) != s

# ---- 個別検証 ----

from rdkit import Chem
_pat = Chem.MolFromSmarts('[CX3](=O)-[OX2]-[#6]')
def hydrolyze_all(smi):
    m = Chem.MolFromSmiles(smi)
    bonds = sorted({m.GetBondBetweenAtoms(a, o).GetIdx() for a, _, o, _ in m.GetSubstructMatches(_pat)})
    if not bonds: return [Chem.MolToSmiles(m)]
    rw = Chem.RWMol(Chem.FragmentOnBonds(m, bonds))
    for at in list(rw.GetAtoms()):
        if at.GetAtomicNum() == 0:
            nb = at.GetNeighbors()[0]; at.SetIsotope(0)
            at.SetAtomicNum(8 if nb.GetSymbol() == 'C' else 1); at.SetNoImplicit(False)
    mm = rw.GetMol(); Chem.SanitizeMol(mm); mm = Chem.RemoveHs(mm)
    return sorted(Chem.MolToSmiles(mm).split('.'))
# 加水分解で立体が保たれることの確認（CIP 記号）
from rdkit.Chem import FindMolChiralCenters
t = 'C[C@H](OC(C)=O)CC'; h = [x for x in hydrolyze_all(t) if 'CC' in x and 'O' in x and '=' not in x][0]
assert FindMolChiralCenters(Chem.MolFromSmiles(t))[0][1] == FindMolChiralCenters(Chem.MolFromSmiles(h))[0][1]
A = 'O=C1OC(C)C(C)OC(=O)C(O)C1O'
assert CalcMolFormula(Chem.MolFromSmiles(A)) == 'C8H12O6'
assert [len(r) for r in Chem.GetSymmSSSR(Chem.MolFromSmiles(A))] == [8]
assert len(stereoisomers('OC(=O)C(O)C(O)C(=O)O')) == 3 and len(stereoisomers('CC(O)C(O)C')) == 3
S = stereoisomers(A); print('A の立体異性体', len(S)); assert len(S) == 10
cls = {}
for s in S:
    parts = hydrolyze_all(s)
    tart = [p for p in parts if 'C(=O)O' in p][0]; diol = [p for p in parts if p != tart][0]
    key = ('Tm' if not is_chiral(tart) else 'T*', 'Dm' if not is_chiral(diol) else 'D*')
    cls.setdefault(key, []).append((s, tart, diol))
    if not is_chiral(s): assert key == ('Tm', 'Dm')
print({k: len(v) for k, v in cls.items()})
assert {k: len(v) for k, v in cls.items()} == {('Tm', 'Dm'): 2, ('Tm', 'D*'): 2, ('T*', 'Dm'): 2, ('T*', 'D*'): 4}
assert sum(not is_chiral(s) for s in S) == 2
a, b = cls[('Tm', 'D*')]; assert mirror(a[0]) == b[0]
# (5) 一方の酒石酸だけ
t1 = cls[('T*', 'Dm')][0][1]
fixed = [s for s in S if t1 in hydrolyze_all(s)]; print('(5)', len(fixed)); assert len(fixed) == 3 and all(is_chiral(s) for s in fixed)
assert all(mirror(x) not in fixed for x in fixed)

print('全検証OK:', ITEM)
