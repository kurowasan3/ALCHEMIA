# 検証スクリプト（グループ g02）  python3 -I この.py
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

def check_all_equations(it, extra=()):
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|…', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]', b): continue
            if not (' + ' in b or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', b)): continue
            e = clean_ce(b)
            r = check_equation(e)
            assert r is None, (b, r)
            print('式OK:', b); n += 1
    for e in extra:
        r = check_equation(e); assert r is None, (e, r); print('式OK(追加):', e); n += 1
    return n

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def check_parts(it):
    k = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, k + 1)], (k, nums)
    print('小問数OK:', k)
ITEM = '/home/claude/chemk/outh/items/ORG-ALC-047/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_parts(it); assert check_all_equations(it) >= 2
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors as D
add = AllChem.ReactionFromSmarts('[CX4;!H0:1][CH1:2]=[O:3].[CX3H1:4]=[O:5]>>[O:3]=[C:2][C:1][C:4][O:5]')
dehyd = AllChem.ReactionFromSmarts('[CX4;!H0:1]([CH1:2]=[O:3])[CX4:4][OH1:5]>>[C:1]([C:2]=[O:3])=[C:4]')
can = lambda m: Chem.MolToSmiles(m, isomericSmiles=False)
def react(*smis):
    ms = [Chem.MolFromSmiles(s) for s in smis]; prods = set()
    for P in ms:
        for Q in ms:
            for ps in add.RunReactants((P, Q)):
                a = ps[0]; Chem.SanitizeMol(a)
                d = dehyd.RunReactants((a,))
                if d:
                    for qs in d:
                        b = qs[0]; Chem.SanitizeMol(b); prods.add(can(b))
                else:
                    prods.add(can(a))
    return prods
f = lambda s: D.CalcMolFormula(Chem.MolFromSmiles(s))
p1 = react('CC(C)C=O'); print(p1, [f(s) for s in p1]); assert [f(s) for s in p1] == ['C8H16O2']
p2 = react('CCC=O', 'CC(C)C=O'); print(p2, sorted(f(s) for s in p2)); assert sorted(f(s) for s in p2) == ['C6H10O', 'C7H12O', 'C7H14O2', 'C8H16O2']
cc = sorted(f(s) for s in p2 if Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts('C=C'))); print('C=C あり', cc); assert cc == ['C6H10O', 'C7H12O']
assert can(Chem.MolFromSmiles('CCC(O)C(C)(C)C=O')) in p2 and can(Chem.MolFromSmiles('CC(C)C=C(C)C=O')) in p2
p3 = react('CCC=O', 'CCCC=O'); fs = sorted(f(s) for s in p3); print(p3, fs); assert len(p3) == 4 and fs.count('C7H12O') == 2
def cleave(smi):
    m = Chem.RWMol(Chem.MolFromSmiles(smi))
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE and b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C':
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx(); break
    m.RemoveBond(i, j)
    for k in (i, j):
        o = m.AddAtom(Chem.Atom(8)); m.AddBond(k, o, Chem.BondType.DOUBLE)
    mm = m.GetMol(); Chem.SanitizeMol(mm)
    return sorted(can(Chem.MolFromSmiles(x)) for x in Chem.MolToSmiles(mm).split('.'))
for X in ['CCCC=O', 'CC(C)C=O']:
    Ys = [s for s in react('CCC=O', X) if f(s) == 'C7H12O']
    print('X', X, 'Y候補', Ys, [cleave(y) for y in Ys])
    assert any(cleave(y) == sorted([can(Chem.MolFromSmiles(X)), can(Chem.MolFromSmiles('CC(=O)C=O'))]) for y in Ys)
    print('X 単独:', [f(s) for s in react(X)])
assert [f(s) for s in react('CCCC=O')] == ['C8H14O'] and [f(s) for s in react('CC(C)C=O')] == ['C8H16O2']
Y = can(Chem.MolFromSmiles('CCCC=C(C)C=O')); assert Y in react('CCC=O', 'CCCC=O') and f(Y) == 'C7H12O'
print('全検証OK:', ITEM)
