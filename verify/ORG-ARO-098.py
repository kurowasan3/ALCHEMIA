# 検証スクリプト（グループ g03）  python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/ORG-ARO-098/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
from rdkit import Chem
from fractions import Fraction as F
from collections import defaultdict
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def run(smi):
    """モデルどおりに生成物を列挙．C_Cl を 14C（同位体 14）にして，生成物ごとに 14C が NH2 の炭素かどうかも記録"""
    m = Chem.MolFromSmiles(smi)
    cl = [a for a in m.GetAtoms() if a.GetSymbol() == 'Cl'][0]
    c = cl.GetNeighbors()[0]
    m = Chem.RWMol(m); m.GetAtomWithIdx(c.GetIdx()).SetIsotope(14)
    orthoH = [n for n in c.GetNeighbors() if n.GetIsAromatic() and n.GetTotalNumHs() > 0]
    out = defaultdict(F)
    for o in orthoH:
        for t in (c.GetIdx(), o.GetIdx()):
            rw = Chem.RWMol(m)
            n = rw.AddAtom(Chem.Atom(7)); rw.AddBond(t, n, Chem.BondType.SINGLE)
            rw.RemoveAtom(cl.GetIdx())
            p = rw.GetMol(); Chem.SanitizeMol(p)
            plain = Chem.Mol(p)
            for a in plain.GetAtoms(): a.SetIsotope(0)
            nC = [a for a in p.GetAtoms() if a.GetSymbol() == 'N'][0].GetNeighbors()[0]
            out[(Chem.MolToSmiles(plain), nC.GetIsotope() == 14)] += F(1, len(orthoH)) * F(1, 2)
    return dict(out)
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
A = {'2,3': can('Cc1cccc(N)c1C'), '3,4': can('Cc1ccc(N)cc1C'), '2,4': can('Cc1ccc(N)c(C)c1'),
     '3,5': can('Cc1cc(C)cc(N)c1'), '2,5': can('Cc1ccc(C)c(N)c1'), '3-Me': can('Cc1cccc(N)c1'), '4-Me': can('Cc1ccc(N)cc1')}
# (1)
r = run('Cc1ccc(Cl)cc1'); print('p-クロロトルエン', r)
assert r == {(A['4-Me'], True): F(1, 2), (A['3-Me'], False): F(1, 2)}
# (2)-(4) 6 種類のクロロキシレン
X = {'2-クロロ-1,3': 'Cc1cccc(C)c1Cl', '2-クロロ-1,4': 'Cc1ccc(C)c(Cl)c1', '4-クロロ-1,2': 'Cc1ccc(Cl)cc1C',
     '3-クロロ-1,2': 'Cc1cccc(Cl)c1C', '4-クロロ-1,3': 'Cc1ccc(Cl)c(C)c1', '5-クロロ-1,3': 'Cc1cc(C)cc(Cl)c1'}
assert len({can(s) for s in X.values()}) == 6 and all(CalcMolFormula(Chem.MolFromSmiles(s)) == 'C8H9Cl' for s in X.values())
def dist(smi):
    d = defaultdict(F)
    for (p, lab), w in run(smi).items(): d[p] += w
    return dict(d)
D = {k: dist(s) for k, s in X.items()}
for k, v in D.items(): print(k, v)
assert D['2-クロロ-1,3'] == {}
assert [k for k, v in D.items() if len(v) == 1] == ['2-クロロ-1,4'] and D['2-クロロ-1,4'] == {A['2,5']: 1}
assert [k for k, v in D.items() if sorted(v.values()) == [F(1, 4), F(3, 4)]] == ['4-クロロ-1,2']
assert D['4-クロロ-1,2'] == {A['3,4']: F(3, 4), A['2,3']: F(1, 4)}
same = [k for k, v in D.items() if v == {A['2,4']: F(1, 2), A['3,5']: F(1, 2)}]
assert same == ['4-クロロ-1,3', '5-クロロ-1,3']
r4 = run(X['4-クロロ-1,3']); r5 = run(X['5-クロロ-1,3'])
print(r4, r5)
assert r4 == {(A['2,4'], True): F(1, 2), (A['3,5'], False): F(1, 2)}
assert r5 == {(A['3,5'], True): F(1, 2), (A['2,4'], False): F(1, 2)}
print('OK')
