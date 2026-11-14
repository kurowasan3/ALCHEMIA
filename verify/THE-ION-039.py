# 検証スクリプト（グループ g09）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out2/items/THE-ION-039/item.json'
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

from sympy import symbols, nsolve, log, sqrt, Rational as R, solve
Ka, Kw, S0, Ksp = 1e-5, 1e-14, 2e-4, 1e-12
import math
def exact_pH(C):
    # 電荷収支: C + h = C*Ka/(h+Ka) + Kw/h を pH について二分法で解く
    f = lambda pH: C + 10**-pH - C*Ka/(10**-pH+Ka) - Kw/10**-pH
    lo, hi = 0.0, 14.0
    for _ in range(200):
        mid = (lo+hi)/2
        if f(lo)*f(mid) <= 0: hi = mid
        else: lo = mid
    return (lo+hi)/2
p1 = exact_pH(1e-2); print('pH(1)', p1); assert round(p1, 1) == 8.5
x = math.sqrt(Kw/Ka*1e-2); assert x/1e-2 < 1e-3 and x > 10*1e-7 and x < S0
h2 = S0*Ka/(1e-2 - S0); print('h2', h2, -math.log10(h2)); assert round(-math.log10(h2), 1) == 6.7
ha7 = 1e-2*1e-7/(1e-7+Ka); print('HA at 7', ha7); assert ha7 < S0 and sig(ha7, 2) == 9.9e-5
a = math.sqrt(Ksp/1e-3); print('A- ppt', a); assert sig(a, 2) == 3.2e-5
need = 1e-2 + 2*(1e-3 - Ksp/1e-4); print('need', need); assert sig(need, 2) == 1.2e-2
p4 = exact_pH(1e-5); print('pH(4)', p4); assert round(p4, 2) == 7.15
print('approx pH(4)', 14 + math.log10(math.sqrt(Kw/Ka*1e-5)))
print('log values', math.log10(2), math.log10(2.04), math.log10(3.16), math.log10(1.41))

print('全検証OK:', ITEM)
