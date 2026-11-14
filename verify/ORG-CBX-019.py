# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out2/items/ORG-CBX-019/item.json'
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
from itertools import product
# (1) 立体異性体の数：P,O,L を アセチル・プロピオニル・ブチリル で代用
acyl = {'P': 'C(=O)C', 'O': 'C(=O)CC', 'L': 'C(=O)CCC'}
def tg(a, b, c):
    return 'C(O%s)C(O%s)CO%s' % (acyl[a], acyl[b], acyl[c])
def count(combo):
    S = set()
    from itertools import permutations
    for p in set(permutations(combo)):
        for s in stereoisomers(tg(*p)):
            S.add(s)
    return S
s1 = count('POL'); print('POL', len(s1)); assert len(s1) == 6
s2 = count('PPO'); ch = [s for s in s2 if is_chiral(s)]
print('PPO', len(s2), 'chiral', len(ch)); assert len(s2) == 3 and len(ch) == 2
assert all(is_chiral(s) for s in s1)
# (2) P-O-L は不斉炭素1つ → 2種
assert len(stereoisomers(tg('P', 'O', 'L'))) == 2
# (3)
x = {'P': R(30, 100), 'O': R(50, 100), 'L': R(20, 100)}
x2 = {'P': R(6, 100), 'O': R(64, 100), 'L': R(30, 100)}
x13 = {k: (3 * x[k] - x2[k]) / 2 for k in x}
print(x13); assert x13 == {'P': R(42, 100), 'O': R(43, 100), 'L': R(15, 100)} and sum(x13.values()) == 1
# (4)
same = sum(v**2 for v in x13.values()); chiral = 1 - same
print('chiral frac', float(chiral)); assert sig(float(chiral) * 100, 2) == 62
# 鏡像の対の確率が等しい（1↔3 入れかえ）
for a, b, c in product('POL', repeat=3):
    assert x13[a] * x2[b] * x13[c] == x13[c] * x2[b] * x13[a]
# (5)
pop = x13['P'] * x2['O'] * x13['P']; print('POP', float(pop)); assert sig(float(pop) * 100, 2) == 11

print('全検証OK:', ITEM)
