# 検証スクリプト（グループ g02）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-ALC-022/item.json'
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
        if ' + ' not in b and not re.search(r'(^|\s)\d+\s*[A-Z(\[]', b): continue
        r = check_equation(b); assert r is None, (b, r); n_eq += 1; print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
from math import log10, floor
def sig(x, n): return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
assert CalcMolFormula(Chem.MolFromSmiles('CC=O')) == CalcMolFormula(Chem.MolFromSmiles('C1CO1')) == 'C2H4O'
assert CalcMolFormula(Chem.MolFromSmiles('OCCO')) == 'C2H6O2'
# OH の数：Z だけが 2 個 → H2 は Z と等モル
oh = lambda s: len(Chem.MolFromSmiles(s).GetSubstructMatches(Chem.MolFromSmarts('[OX2H]')))
assert (oh('CC=O'), oh('C1CO1'), oh('OCCO')) == (0, 0, 2)
print('Z 0.10 mol → H2', 0.10 * oh('OCCO') / 2, 'mol')
print('OK')
