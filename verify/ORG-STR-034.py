# 検証スクリプト ORG-STR-034（グループ g01）  python3 -I この.py
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

def check_item(path):
    it = json.load(open(path, encoding='utf-8'))
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
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    return it

it = check_item('/home/claude/chem5/out/items/ORG-STR-034/item.json')

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
R = {'a': 'C', 'b': 'CC', 'c': 'c1ccccc1', 'd': 'c1ccccc1C'}
G = {'p': '', 'q': 'O', 'r': 'C=O', 's': 'C(C)=O', 't': 'OC', 'u': 'C(=O)O', 'v': 'C'}
d = {}
for r, rs in R.items():
    for g, gs in G.items():
        c = Chem.MolToSmiles(Chem.MolFromSmiles(rs + gs)); d.setdefault(c, []).append(r + g)
dups = {c: sorted(v) for c, v in d.items() if len(v) > 1}
print(dups)
assert dups == {'CC': ['av', 'bp'], 'Cc1ccccc1': ['cv', 'dp']} and len(d) == 26
F = {k: CalcMolFormula(Chem.MolFromSmiles(R[k[0]] + G[k[1]])) for k in ['br', 'as', 'cs', 'dq', 'ct']}
print(F)
assert F == {'br': 'C3H6O', 'as': 'C3H6O', 'cs': 'C8H8O', 'dq': 'C7H8O', 'ct': 'C7H8O'}

print('ALL OK')
