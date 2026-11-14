# 検証スクリプト ORG-CBX-042（グループ g09）  python3 -I この.py
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

it = check_item('/home/claude/chem6/out/items/ORG-CBX-042/item.json')
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as F
S = {'P':'CC(=O)OCC(O)C', 'Q':'CC(OC(C)=O)CO', "P'":'CC(=O)OCC(C)=O', "Q'":'CC(OC(C)=O)C=O',
     'Pdi':'CC(=O)OCC(C)OC(C)=O', 'Qdi':'CC(OC(C)=O)COC(C)=O'}
for k, s in S.items():
    m = Chem.MolFromSmiles(s)
    cc = Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
    print(k, F(m), '不斉炭素', len(cc))
chir = {k: len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False)) for k, s in S.items()}
assert [k for k in ('P','Q',"P'","Q'") if chir[k]] == ['P','Q',"Q'"]
assert F(Chem.MolFromSmiles(S['P'])) == 'C5H10O3' == F(Chem.MolFromSmiles(S['Q']))
assert Chem.CanonSmiles(S['Pdi']) == Chem.CanonSmiles(S['Qdi']) and F(Chem.MolFromSmiles(S['Pdi'])) == 'C7H12O4'
# 銀鏡: P' はケトン, Q' はアルデヒド
ald = Chem.MolFromSmarts('[CX3H1](=O)[#6]')
assert not Chem.MolFromSmiles(S["P'"]).HasSubstructMatch(ald) and Chem.MolFromSmiles(S["Q'"]).HasSubstructMatch(ald)
M = 5*12 + 10*1.0 + 3*16; n = 1.18 / M; V = n / 0.500 * 1000
print('M =', M, 'n =', n, 'V =', V); assert abs(V - 20.0) < 1e-9
print('OK')
