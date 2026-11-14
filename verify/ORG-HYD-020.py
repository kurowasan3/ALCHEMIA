# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out2/items/ORG-HYD-020/item.json'
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

from rdkit.Chem import AllChem
c4h6 = ['C#CCC', 'CC#CC', 'C=C=CC', 'C=CC=C']
assert all(CalcMolFormula(Chem.MolFromSmiles(s)) == 'C4H6' for s in c4h6)
assert [s for s in c4h6 if Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts('C=C-C=C'))] == ['C=CC=C']
da = AllChem.ReactionFromSmarts('[C:1]=[C:2]-[C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1')
iso = Chem.MolFromSmiles('C=C(C)C=C')
P = set()
for r in da.RunReactants((iso, iso)):
    p = r[0]; Chem.SanitizeMol(p); P.add(Chem.MolToSmiles(p))
print(sorted(P)); assert len(P) == 4 and all(CalcMolFormula(Chem.MolFromSmiles(p)) == 'C10H16' for p in P)
assert sum(len(stereoisomers(p)) for p in P) == 8
oz = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>([C:1]=O.[C:2]=O)')
def ozall(s):
    todo, done = [s], []
    while todo:
        f = todo.pop(); ps = oz.RunReactants((Chem.MolFromSmiles(f),))
        if not ps: done.append(Chem.MolToSmiles(Chem.MolFromSmiles(f))); continue
        m = ps[0][0]; Chem.SanitizeMol(m); todo += Chem.MolToSmiles(m).split('.')
    return sorted(done)
cho = Chem.MolFromSmarts('[CX3H1](=O)[#6]')
def hydro(s):
    m = Chem.RWMol(Chem.MolFromSmiles(s))
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE: b.SetBondType(Chem.BondType.SINGLE)
    mm = m.GetMol(); Chem.SanitizeMol(mm); return Chem.MolToSmiles(mm)
cand = []
for p in sorted(P):
    prods = ozall(p); main = [x for x in prods if x != 'C=O']
    assert len(main) == 1 and prods.count('C=O') == 1
    n = main[0]; ncho = len(Chem.MolFromSmiles(n).GetSubstructMatches(cho))
    M = hydro(p); Ms = stereoisomers(M); allinact = all(not is_chiral(x) for x in Ms)
    print(p, n, CalcMolFormula(Chem.MolFromSmiles(n)), 'CHO', ncho, 'M', M, len(Ms), allinact)
    if CalcMolFormula(Chem.MolFromSmiles(n)) == 'C9H14O3' and ncho == 1 and allinact: cand.append((p, n, len(Ms)))
print(cand); assert len(cand) == 1 and cand[0][0] == Chem.MolToSmiles(Chem.MolFromSmiles('C=C(C)C1CC=C(C)CC1')) and cand[0][2] == 2
assert Chem.MolFromSmiles(cand[0][1]).HasSubstructMatch(Chem.MolFromSmarts('[CH3]C(=O)[#6]'))
assert all(is_chiral(x) for x in stereoisomers(cand[0][1]))

print('全検証OK:', ITEM)
