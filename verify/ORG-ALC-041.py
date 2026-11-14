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
ITEM = '/home/claude/chemk/out/items/ORG-ALC-041/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_parts(it); assert check_all_equations(it) >= 3
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1]([H])[C:2]C=O')
def hf(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    for a, b in [(0, 1), (1, 0)]:
        pass
    r1 = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1][C:2]C=O')
    for ps in r1.RunReactants((m,)):
        p = ps[0]; Chem.SanitizeMol(p)
        out.add(Chem.MolToSmiles(Chem.RemoveHs(p), isomericSmiles=False))
    return out
alk = {'1-ブテン': 'C=CCC', 'シス-2-ブテン': 'C/C=C\\C', 'トランス-2-ブテン': 'C/C=C/C', '2-メチルプロペン': 'C=C(C)C'}
res = {k: hf(v) for k, v in alk.items()}
for k, v in res.items(): print(k, v)
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=False)
assert res['1-ブテン'] == {can('CCCCC=O'), can('CCC(C)C=O')}
assert res['シス-2-ブテン'] == res['トランス-2-ブテン'] == {can('CCC(C)C=O')}
assert res['2-メチルプロペン'] == {can('CC(C)CC=O'), can('CC(C)(C)C=O')}
single = [k for k, v in res.items() if len(v) == 1]; assert single == ['シス-2-ブテン', 'トランス-2-ブテン']
src = [k for k, v in res.items() if can('CC(C)CC=O') in v]; assert src == ['2-メチルプロペン']
for s in ['CCCCC=O', 'CCC(C)C=O']:
    print(s, Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True))
assert len(Chem.FindMolChiralCenters(Chem.MolFromSmiles('CCC(C)C=O'), includeUnassigned=True)) == 1
assert len(Chem.FindMolChiralCenters(Chem.MolFromSmiles('CCCCC=O'), includeUnassigned=True)) == 0
mw = lambda f: sum({'C': 12, 'H': 1.0, 'O': 16}[e]*n for e, n in f)
anh = Chem.MolFromSmiles('O=C1OC(=O)c2ccccc12'); est = Chem.MolFromSmiles('CCCCCOC(=O)c1ccccc1C(=O)OCCCCC')
print(rdMolDescriptors.CalcMolFormula(anh), rdMolDescriptors.CalcMolFormula(est))
assert rdMolDescriptors.CalcMolFormula(anh) == 'C8H4O3' and rdMolDescriptors.CalcMolFormula(est) == 'C18H26O4'
Manh = 8*12 + 4*1.0 + 3*16; Mest = 18*12 + 26*1.0 + 4*16
n = 29.6/Manh; print(Manh, Mest, n, Mest*n); assert sig(Mest*n, 3) == 61.2
but = 2*n/0.80*56; print('1-butene g', but); assert sig(but, 3) == 28.0
print('全検証OK:', ITEM)
