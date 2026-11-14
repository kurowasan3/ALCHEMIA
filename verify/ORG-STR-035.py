# 検証スクリプト ORG-STR-035（グループ g01）  python3 -I この.py
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

it = check_item('/home/claude/chem5/out/items/ORG-STR-035/item.json')

from rdkit import Chem
S = {'ア': 'C(=O)[OH]', 'イ': 'c[OH]', 'ウ': '[CX4][OH]', 'エ': '[#6]C(=O)O[#6]', 'オ': '[N+](=O)[O-]', 'カ': '[NX3;H2][CX4]C(=O)[OH]', 'キ': 'c'}
M = {1: 'OC(=O)c1ccccc1O', 2: 'Oc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]', 3: 'CC(O)C(=O)O', 4: 'CC(N)C(=O)O',
     5: 'CC(=O)Oc1ccccc1C(=O)O', 6: 'OCc1ccccc1', 7: 'Cc1ccccc1O', 8: 'COC(=O)c1ccccc1O'}
ANS = {1: 'アイキ', 2: 'イオキ', 3: 'アウ', 4: 'アカ', 5: 'アエキ', 6: 'ウキ', 7: 'イキ', 8: 'イエキ'}
for i, s in M.items():
    m = Chem.MolFromSmiles(s)
    got = ''.join(k for k, p in S.items() if m.HasSubstructMatch(Chem.MolFromSmarts(p)))
    print(i, got); assert got == ANS[i], (i, got)
fecl3 = [i for i in (1, 5, 6, 7, 8) if Chem.MolFromSmiles(M[i]).HasSubstructMatch(Chem.MolFromSmarts(S['イ']))]
assert fecl3 == [1, 7, 8]

print('ALL OK')
