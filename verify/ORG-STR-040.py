# 検証スクリプト（グループ g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-STR-040/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)

# ---- 計算
from sympy import Rational as R, symbols, solve
H, C, O, Br = R(1), R(12), R(16), R(80)
nBr2 = R(40)*R(1,10)/(2*Br)
MX = R(140,100)/nBr2
assert MX == 56
n = MX/14; assert n == 4
a, b = symbols('a b')
sol = solve([4*a+7*b-R(792,10)/44, 10*a+8*b-2*R(234,10)/18], [a, b])
print('A, トルエン (mmol) =', sol)
assert sol[a] == R(1,10) and sol[b] == R(2,10)
MA = 4*C+10*H+O; MT = 7*C+8*H
mass = MA*sol[a] + MT*sol[b]
assert mass == R(258,10), mass
pct = float(MA*sol[a]/mass*100); print('A の質量% =', pct)
assert round(pct) == 29

# ---- 構造：RDKit で C4H8 アルケンの水付加（マルコフニコフ）と C4H10O アルコールの脱水を網羅
from rdkit import Chem
from rdkit.Chem import AllChem
alkenes = ['C=CCC', 'CC=CC', 'C=C(C)C']
def markov(smi):
    m = Chem.MolFromSmiles(smi)
    for bd in m.GetBonds():
        if bd.GetBondTypeAsDouble() == 2:
            a1, a2 = bd.GetBeginAtom(), bd.GetEndAtom()
            # OH は H の少ない炭素へ
            tgt = a1 if a1.GetTotalNumHs() <= a2.GetTotalNumHs() else a2
            if a1.GetTotalNumHs() == a2.GetTotalNumHs():
                return None
            rw = Chem.RWMol(m); rw.GetBondWithIdx(bd.GetIdx()).SetBondType(Chem.BondType.SINGLE)
            o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(tgt.GetIdx(), o, Chem.BondType.SINGLE)
            mm = rw.GetMol(); Chem.SanitizeMol(mm); return Chem.MolToSmiles(mm)
def degree(smi):
    m = Chem.MolFromSmiles(smi)
    for a in m.GetAtoms():
        if a.GetSymbol() == 'O':
            c = a.GetNeighbors()[0]; return sum(1 for x in c.GetNeighbors() if x.GetSymbol() == 'C')
tert = [s for s in alkenes if markov(s) and degree(markov(s)) == 3]
print('水付加で第三級アルコールを与える C4H8:', tert)
assert tert == ['C=C(C)C']
dehyd = AllChem.ReactionFromSmarts('[C;!H0:1][C:2][OX2H1]>>[C:1]=[C:2]')
def dehydrate(smi):
    m = Chem.MolFromSmiles(smi)
    ps = set()
    for p in dehyd.RunReactants((m,)):
        q = p[0]; Chem.SanitizeMol(q); ps.add(Chem.MolToSmiles(q))
    return ps
alcs = ['OCCCC', 'CC(O)CC', 'OCC(C)C', 'CC(C)(C)O']
ok = [s for s in alcs if dehydrate(s) == {Chem.CanonSmiles('C=C(C)C')} and degree(s) == 1]
print('脱水で X だけを与える第一級アルコール:', ok)
assert ok == ['OCC(C)C']
print('OK')
