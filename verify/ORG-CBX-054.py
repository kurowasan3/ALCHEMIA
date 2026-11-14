# 検証スクリプト ORG-CBX-054  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/out/items/ORG-CBX-054/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    if n_items: assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()


from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
A = Chem.MolFromSmiles('CCOC(=O)/C=C\\C(=O)OC(C)CC')
B = Chem.MolFromSmiles('CCOC(=O)/C=C/C(=O)OC(C)(C)C')
for m in (A, B): assert CalcMolFormula(m) == 'C10H16O4'
nc = lambda m: len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
assert nc(A) == 1 and nc(B) == 0
assert CalcMolFormula(Chem.MolFromSmiles('OC(=O)/C=C\\C(=O)O')) == 'C4H4O4'
for s in ('CC(O)CC', 'CC(C)(C)O'): assert CalcMolFormula(Chem.MolFromSmiles(s)) == 'C4H10O'
C = 35.2*12/44; H = 18.0*2.0/18; O = 14.8 - C - H
r = [C/12, H/1.0, O/16]; m = min(r); print([x/m for x in r])
assert [round(x/m) for x in r] == [4, 10, 1]
print('OK')
