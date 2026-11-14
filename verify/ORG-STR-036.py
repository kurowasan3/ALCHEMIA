# 検証スクリプト ORG-STR-036（グループ g01）  python3 -I この.py
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

it = check_item('/home/claude/chem5/out/items/ORG-STR-036/item.json')

from rdkit import Chem
def P(smi):
    m = Chem.MolFromSmiles(smi)
    has = lambda q: m.HasSubstructMatch(Chem.MolFromSmarts(q))
    a = len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)) > 0
    b = has('[CX3]=[CX3]') and not has('c')  or has('[CH2]=[CH]c')
    c = len(m.GetSubstructMatches(Chem.MolFromSmarts('C(=O)[OH]'))) == 2
    d = has('[CH3][CX3](=O)[#6,#1]') or has('[CH3][CX3H1]=O') or has('[CH3][CX4;H1,H2][OH]')
    e = has('[CX3H1](=O)') or has('[CH2]=O')
    f = has('[#6]C(=O)O[#6]')
    g = has('[NX3;H2][CX4]C(=O)[OH]')
    return dict(a=a, b=b, c=c, d=d, e=e, f=f, g=g)
groups = [
 ({'乳酸': 'CC(O)C(=O)O', 'アラニン': 'CC(N)C(=O)O', '2-ブタノール': 'CCC(C)O', '2-クロロブタン': 'CCC(C)Cl'}, ('グリシン', 'NCC(=O)O'), 'a'),
 ({'スチレン': 'C=Cc1ccccc1', '酢酸ビニル': 'C=COC(C)=O', 'アクリル酸メチル': 'C=CC(=O)OC', 'プロペン': 'C=CC'}, ('エチレングリコール', 'OCCO'), 'b'),
 ({'シュウ酸': 'OC(=O)C(=O)O', 'マレイン酸': 'OC(=O)/C=C\\C(=O)O', 'フタル酸': 'OC(=O)c1ccccc1C(=O)O', 'アジピン酸': 'OC(=O)CCCCC(=O)O'}, ('乳酸', 'CC(O)C(=O)O'), 'c'),
 ({'エタノール': 'CCO', 'アセトン': 'CC(C)=O', '2-プロパノール': 'CC(C)O', 'アセトアルデヒド': 'CC=O'}, ('1-プロパノール', 'CCCO'), 'd'),
 ({'ホルムアルデヒド': 'C=O', 'ギ酸': 'O=CO', 'グルコース': 'OCC(O)C(O)C(O)C(O)C=O', 'アセトアルデヒド': 'CC=O'}, ('アセトン', 'CC(C)=O'), 'e'),
]
for four, (oname, osmi), key in groups:
    props = [P(s) for s in four.values()]; po = P(osmi)
    common = [k for k in 'abcdefg' if all(p[k] for p in props) and not po[k]]
    print(key, common, oname)
    assert common == [key], (key, common)

print('ALL OK')
