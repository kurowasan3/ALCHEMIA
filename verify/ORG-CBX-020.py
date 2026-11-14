# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out2/items/ORG-CBX-020/item.json'
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

from sympy import Rational as R, symbols, solve
from rdkit.Chem import AllChem
U = {'U1': 'CCCCCCCCC=CCCCCCCCC(=O)O', 'U2': 'CCCCCC=CCC=CCCCCCCCC(=O)O', 'U3': 'CCC=CCC=CCC=CCCCCCCCC(=O)O'}
fml = {'U1': 'C18H34O2', 'U2': 'C18H32O2', 'U3': 'C18H30O2'}
rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>([C:1](=O)O.[C:2](=O)O)')
def cleave(smi):
    frags = [smi]
    while True:
        new = []; changed = False
        for f in frags:
            ps = rxn.RunReactants((Chem.MolFromSmiles(f),))
            if ps:
                p = ps[0][0]; Chem.SanitizeMol(p); new += Chem.MolToSmiles(p).split('.'); changed = True
            else: new.append(f)
        frags = new
        if not changed: return sorted(Chem.MolToSmiles(Chem.MolFromSmiles(x)) for x in frags)
c = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
aze, mal, non, hex_, pro = c('OC(=O)CCCCCCCC(=O)O'), c('OC(=O)CC(=O)O'), c('CCCCCCCCC(=O)O'), c('CCCCCC(=O)O'), c('CCC(=O)O')
exp = {'U1': sorted([non, aze]), 'U2': sorted([hex_, mal, aze]), 'U3': sorted([pro, mal, mal, aze])}
for k in U:
    assert CalcMolFormula(Chem.MolFromSmiles(U[k])) == fml[k]
    print(k, cleave(U[k])); assert cleave(U[k]) == exp[k]
x = {'S': R(10, 100), 'U1': R(25, 100), 'U2': R(45, 100), 'U3': R(20, 100)}
nb = x['U1'] + 2*x['U2'] + 3*x['U3']; a = x['U1'] + x['U2'] + x['U3']; m = x['U2'] + 2*x['U3']; p = x['U3']
print('ratio', 3*a, 3*m, 3*p); assert (3*a, 3*m, 3*p) == (R(270, 100), R(255, 100), R(60, 100))
M_tristearin = 57*12 + 110*R(1) + 6*16; assert M_tristearin == 890
IV = 3*nb*254/(M_tristearin - 6*nb)*100; print('IV', float(IV)); assert round(float(IV)) == 152
n = symbols('n'); nsol = solve(3*n*254/(890 - 6*n)*100 - 152, n)[0]; print('nbar', float(nsol)); assert sig(float(nsol), 3) == 1.75
s, x1, x2, x3 = symbols('s x1 x2 x3')
sol = solve([s+x1+x2+x3-1, x1+x2+x3-R(90, 100), x2+2*x3-R(85, 100), x3-R(20, 100)], [s, x1, x2, x3])
print(sol); assert sol == {s: R(1, 10), x1: R(1, 4), x2: R(9, 20), x3: R(1, 5)}
assert a + m == nb
t = R(5, 100); y = {'S': x['S'], 'U1': x['U1']+t, 'U2': x['U2']-2*t, 'U3': x['U3']+t}
assert y['U1']+y['U2']+y['U3'] == a and y['U2']+2*y['U3'] == m and y['U1']+2*y['U2']+3*y['U3'] == nb and y['U3'] != p
print('別解', y)

print('全検証OK:', ITEM)
