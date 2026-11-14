# 検証スクリプト ORG-CBX-041（グループ g09）  python3 -I この.py
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
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
    assert np_ == na or (np_ == 0 and na == 0), ('小問の数', np_, na)
    print('小問の数', np_, na)
    return it

it = check_item('/home/claude/chem6/out/items/ORG-CBX-041/item.json')
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
smi = {'酢酸':'CC(=O)O','エタノール':'CCO','酢酸エチル':'CCOC(C)=O','フェノール':'Oc1ccccc1','無水酢酸':'CC(=O)OC(C)=O',
       '酢酸フェニル':'CC(=O)Oc1ccccc1','アセチレン':'C#C','酢酸ビニル':'C=COC(C)=O','サリチル酸':'O=C(O)c1ccccc1O',
       'メタノール':'CO','サリチル酸メチル':'COC(=O)c1ccccc1O','アセチルサリチル酸':'CC(=O)Oc1ccccc1C(=O)O','水':'O'}
from collections import Counter
def atoms(n):
    m = Chem.AddHs(Chem.MolFromSmiles(smi[n])); return Counter(a.GetSymbol() for a in m.GetAtoms())
def bal(L, R):
    a = sum((atoms(x) for x in L), Counter()); b = sum((atoms(x) for x in R), Counter()); assert a == b, (L, R, a, b)
bal(['酢酸','エタノール'], ['酢酸エチル','水'])
bal(['無水酢酸','フェノール'], ['酢酸フェニル','酢酸'])
bal(['酢酸','アセチレン'], ['酢酸ビニル'])
bal(['サリチル酸','メタノール'], ['サリチル酸メチル','水'])
bal(['サリチル酸','無水酢酸'], ['アセチルサリチル酸','酢酸'])
# 生成物の官能基: アセチルサリチル酸はCOOHが残りフェノール性OHなし / サリチル酸メチルはフェノール性OHが残る
asp = Chem.MolFromSmiles(smi['アセチルサリチル酸']); ms = Chem.MolFromSmiles(smi['サリチル酸メチル'])
phenol = Chem.MolFromSmarts('c[OX2H]'); cooh = Chem.MolFromSmarts('C(=O)[OX2H]')
assert not asp.HasSubstructMatch(phenol) and asp.HasSubstructMatch(cooh)
assert ms.HasSubstructMatch(phenol) and not ms.HasSubstructMatch(cooh)
print('原子の収支・官能基 OK')
