# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out2/items/ORG-CBX-021/item.json'
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

M = {'H': 1.0, 'C': 12, 'O': 16}
benzoate = 8*12 + 8*1.0 + 2*16; assert benzoate == 136
assert benzoate + 2 == 138 and 2*1.0 + 18 == 20 and benzoate + 4 == 140
# 立体：アシル-O 開裂で不斉炭素の配置が保たれる
from rdkit.Chem import FindMolChiralCenters
_pat = Chem.MolFromSmarts('[CX3](=O)-[OX2]-[#6]')
def acyl_cleave(smi):
    m = Chem.MolFromSmiles(smi)
    bonds = sorted({m.GetBondBetweenAtoms(a, o).GetIdx() for a, _, o, _ in m.GetSubstructMatches(_pat)})
    rw = Chem.RWMol(Chem.FragmentOnBonds(m, bonds))
    for at in list(rw.GetAtoms()):
        if at.GetAtomicNum() == 0:
            nb = at.GetNeighbors()[0]; at.SetIsotope(0); at.SetAtomicNum(8 if nb.GetSymbol() == 'C' else 1); at.SetNoImplicit(False)
    mm = rw.GetMol(); Chem.SanitizeMol(mm); return Chem.MolToSmiles(Chem.RemoveHs(mm)).split('.')
for alc in ['CC[C@H](C)O', 'CC[C@@H](C)O']:
    ester = can(Chem.MolFromSmiles(alc.replace('O', 'OC(C)=O') if alc.endswith('O') else alc))
    back = [x for x in acyl_cleave(ester) if 'CC' in x and '=O' not in x][0]
    c1 = FindMolChiralCenters(Chem.MolFromSmiles(alc)); c2 = FindMolChiralCenters(Chem.MolFromSmiles(back))
    print(alc, ester, back, c1, c2); assert c1 == c2
assert len(stereoisomers('CCCC(C)(O)CC')) == 2 and len(stereoisomers('CCC(C)O')) == 2

print('全検証OK:', ITEM)
