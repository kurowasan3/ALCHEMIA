# 検証スクリプト（グループ g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-ALC-024/item.json'
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

# ---- RDKit：C5H12O アルコールに「脱水（ザイツェフ）→ 水付加（マルコフニコフ）」を当てはめる
from rdkit import Chem
from rdkit.Chem import AllChem
dehyd = AllChem.ReactionFromSmarts('[C;!H0:1][C:2][OX2H1]>>[C:1]=[C:2]')
def can(s): return Chem.CanonSmiles(s)
def degree(m, oidx=None):
    for a in m.GetAtoms():
        if a.GetSymbol() == 'O':
            c = a.GetNeighbors()[0]; return sum(1 for x in c.GetNeighbors() if x.GetSymbol() == 'C')
def dbH(smi):
    m = Chem.MolFromSmiles(smi)
    for b in m.GetBonds():
        if b.GetBondTypeAsDouble() == 2:
            return b.GetBeginAtom().GetTotalNumHs() + b.GetEndAtom().GetTotalNumHs()
def zaitsev(smi):
    m = Chem.MolFromSmiles(smi); ps = set()
    for p in dehyd.RunReactants((m,)):
        q = p[0]; Chem.SanitizeMol(q); ps.add(Chem.MolToSmiles(Chem.RemoveHs(q)))
    if not ps: return set(), None
    # シス・トランスは区別しない（立体なしの SMILES）
    mn = min(dbH(s) for s in ps)
    return ps, {s for s in ps if dbH(s) == mn}
def markov(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    for bd in m.GetBonds():
        if bd.GetBondTypeAsDouble() == 2:
            a1, a2 = bd.GetBeginAtom(), bd.GetEndAtom()
            h1, h2 = a1.GetTotalNumHs(), a2.GetTotalNumHs()
            tg = [a1] if h1 < h2 else [a2] if h2 < h1 else [a1, a2]
            for t in tg:
                rw = Chem.RWMol(m); rw.GetBondWithIdx(bd.GetIdx()).SetBondType(Chem.BondType.SINGLE)
                o = rw.AddAtom(Chem.Atom(8)); rw.AddBond(t.GetIdx(), o, Chem.BondType.SINGLE)
                mm = rw.GetMol(); Chem.SanitizeMol(mm); out.add(Chem.MolToSmiles(mm))
    return out
alcs = {'1-pentanol':'OCCCCC','2-pentanol':'CC(O)CCC','3-pentanol':'CCC(O)CC','2-methyl-1-butanol':'OCC(C)CC',
        '3-methyl-1-butanol':'OCCC(C)C','3-methyl-2-butanol':'CC(O)C(C)C','2-methyl-2-butanol':'CC(C)(O)CC','2,2-dimethyl-1-propanol':'OCC(C)(C)C'}
toTert = []
for name, s in alcs.items():
    allp, main = zaitsev(s)
    if not main:
        print(name, '脱水できない'); continue
    prods = set().union(*[markov(x) for x in main])
    degs = {degree(Chem.MolFromSmiles(p)) for p in prods}
    print(name, '脱水:', allp, '主:', main, '→ 付加:', prods, '級数:', degs)
    if degree(Chem.MolFromSmiles(s)) != 3 and degs == {3}: toTert.append(name)
print('第三級に変わる:', toTert)
assert sorted(toTert) == ['2-methyl-1-butanol', '3-methyl-2-butanol']
allp, main = zaitsev(alcs['3-methyl-2-butanol'])
assert allp == {can('CC=C(C)C'), can('C=CC(C)C')} and main == {can('CC=C(C)C')}
assert markov('CC=C(C)C') == {can('CCC(C)(C)O')}
_, m1 = zaitsev(alcs['1-pentanol']); assert m1 == {can('C=CCCC')}
assert markov('C=CCCC') == {can('CC(O)CCC')}
print('OK')
