# 検証スクリプト（グループ g03）  python3 -I この.py

import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

def check_all_ce(it):
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if should_check(b):
                r = check_equation(clean_ce(b))
                assert r is None, (b, r)
                print('式OK:', b); n += 1
    print('検査した反応式:', n)

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import combinations
it = json.load(open('/home/claude/chem5/outh/items/ORG-CBX-037/item.json', encoding='utf-8'))
check_all_ce(it)
def can(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s))
cons = set()
for sk in ('CCCC', 'CC(C)C'):
    base = Chem.MolFromSmiles(sk)
    hs = [(a.GetIdx(), k) for a in base.GetAtoms() for k in range(a.GetTotalNumHs())]
    for (i, _), (j, _) in combinations(hs, 2):
        rw = Chem.RWMol(base)
        for t in (i, j):
            c = rw.AddAtom(Chem.Atom(6)); o1 = rw.AddAtom(Chem.Atom(8)); o2 = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(t, c, Chem.BondType.SINGLE); rw.AddBond(c, o1, Chem.BondType.DOUBLE); rw.AddBond(c, o2, Chem.BondType.SINGLE)
        m = rw.GetMol(); Chem.SanitizeMol(m); cons.add(Chem.MolToSmiles(m))
cons = sorted(cons)
assert len(cons) == 9 and all(CalcMolFormula(Chem.MolFromSmiles(s)) == 'C6H10O4' for s in cons)
st = sorted({Chem.MolToSmiles(x) for s in cons for x in EnumerateStereoisomers(Chem.MolFromSmiles(s))})
assert len(st) == 13
print('構造異性体 9，立体異性体 13')
def nC(s): return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
def mirror(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('@@', '!').replace('@', '@@').replace('!', '@')))
def cc_of(m): return [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'C' and sum(1 for n in a.GetNeighbors() if n.GetSymbol() == 'O') == 2]
def heat(s):
    m = Chem.MolFromSmiles(s); c1, c2 = cc_of(m); path = Chem.GetShortestPath(m, c1, c2)
    if len(path) == 3:
        rw = Chem.RWMol(m)
        for n in sorted([n.GetIdx() for n in m.GetAtomWithIdx(c2).GetNeighbors() if n.GetSymbol() == 'O'], reverse=True): rw.RemoveAtom(n)
        rw.RemoveAtom(c2); p = rw.GetMol(); Chem.SanitizeMol(p); return ('CO2', Chem.MolToSmiles(p))
    ring = len(path) + 1
    if ring in (5, 6):
        rw = Chem.RWMol(m)
        oh = lambda c: [n.GetIdx() for n in m.GetAtomWithIdx(c).GetNeighbors() if n.GetSymbol() == 'O' and m.GetBondBetweenAtoms(c, n.GetIdx()).GetBondType() == Chem.BondType.SINGLE][0]
        o1, o2 = oh(c1), oh(c2)
        rw.AddBond(o1, c2, Chem.BondType.SINGLE); rw.RemoveAtom(o2); p = rw.GetMol(); Chem.SanitizeMol(p)
        return ('無水物', ring, Chem.MolToSmiles(p))
    return ('なし', ring)
H = {s: heat(s) for s in cons}
from collections import Counter
cls = Counter((h[0], h[1]) if h[0] != 'CO2' else ('CO2',) for h in H.values())
print(cls)
assert cls[('無水物', 5)] == 3 and cls[('無水物', 6)] == 2 and cls[('CO2',)] == 3 and cls[('なし', 7)] == 1
assert [s for s in cons if H[s][0] == 'なし'] == [can('OC(=O)CCCCC(=O)O')]
# (3) P：不斉炭素 2 個・鏡像と同一（メソ）
P = [s for s in st if nC(s) == 2 and mirror(s) == s]
assert len(P) == 1 and can(Chem.MolToSmiles(Chem.MolFromSmiles(P[0]), isomericSmiles=False)) == can('CC(C(=O)O)C(C)C(=O)O')
Pp = [s for s in st if nC(s) == 2 and mirror(s) != s]; assert len(Pp) == 2
aP = heat(P[0])[2]; aPp = [heat(s)[2] for s in Pp]
print('P の無水物', aP, 'P′ の無水物', aPp)
assert mirror(aP) == aP and all(mirror(a) != a for a in aPp) and mirror(aPp[0]) == aPp[1]
# 図の立体式
assert can('C[C@H](C(=O)O)[C@@H](C)C(=O)O') == P[0]
assert can('C[C@@H]1C(=O)OC(=O)[C@@H]1C') == aP and can('C[C@H]1C(=O)OC(=O)[C@@H]1C') in aPp
# (4) Q：不斉炭素なし・脱炭酸生成物に不斉炭素
Q = [s for s in cons if nC(s) == 0 and H[s][0] == 'CO2' and nC(H[s][1]) > 0]
assert Q == [can('CCC(C)(C(=O)O)C(=O)O')] and H[Q[0]][1] == can('CCC(C)C(=O)O')
# Q の 2 個の COOH は対称的に等価（正規化ランクが同じ）→ ラセミ体
mQ = Chem.MolFromSmiles(Q[0]); ranks = list(Chem.CanonicalRankAtoms(mQ, breakTies=False))
c1, c2 = cc_of(mQ); assert ranks[c1] == ranks[c2]
# (5) S：不斉炭素 1 個・六員環
S = [s for s in cons if nC(s) == 1 and H[s][:2] == ('無水物', 6)]
assert S == [can('OC(=O)C(C)CCC(=O)O')]
assert check_equation('C6H10O4 -> C6H8O3 + H2O') is None
print('すべて一致')
