# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out2/items/ORG-STR-014/item.json'
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

from sympy import Rational as R
E = 'CC(C(=O)OC(C)c1ccccc1)c1ccccc1'
assert CalcMolFormula(Chem.MolFromSmiles(E)) == 'C17H18O2'
assert 17*12 + 18*1.0 + 2*16 == 254
S = stereoisomers(E); assert len(S) == 4
_pat = Chem.MolFromSmarts('[CX3](=O)-[OX2]-[#6]')
def parts(smi):
    m = Chem.MolFromSmiles(smi)
    bonds = sorted({m.GetBondBetweenAtoms(a, o).GetIdx() for a, _, o, _ in m.GetSubstructMatches(_pat)})
    rw = Chem.RWMol(Chem.FragmentOnBonds(m, bonds))
    for at in list(rw.GetAtoms()):
        if at.GetAtomicNum() == 0:
            nb = at.GetNeighbors()[0]; at.SetIsotope(0); at.SetAtomicNum(8 if nb.GetSymbol() == 'C' else 1); at.SetNoImplicit(False)
    mm = rw.GetMol(); Chem.SanitizeMol(mm); fr = Chem.MolToSmiles(Chem.RemoveHs(mm)).split('.')
    acid = [f for f in fr if 'C(=O)O' in f][0]; alc = [f for f in fr if f != acid][0]
    return acid, alc
for s in S:
    a, b = parts(s); ma, mb = parts(mirror(s))
    assert ma == mirror(a) and mb == mirror(b)   # 鏡像は両方の配置を逆にしたもの
    assert mirror(s) in S and mirror(s) != s
# A(+) 固定で生じる 2 種は鏡像でない
a0 = parts(S[0])[0]
fixed = [s for s in S if parts(s)[0] == a0]; assert len(fixed) == 2 and mirror(fixed[0]) != fixed[1]
# (4)
Bp = Bm = R(1, 20); Ap, Am = R(9, 100), R(1, 100)
ApBm, ApBp, AmBp, AmBm = Bm*R(9, 10), Bp*R(9, 10), Bp*R(1, 10), Bm*R(1, 10)
frac = ApBm + AmBp; mass = frac*254; print('mass', float(mass)); assert sig(float(mass), 3) == 12.7
fracBm = ApBm/frac; rot = (ApBm - AmBp)/frac
print(float(fracBm), float(rot)); assert fracBm == R(9, 10) and rot == R(4, 5)

print('全検証OK:', ITEM)
