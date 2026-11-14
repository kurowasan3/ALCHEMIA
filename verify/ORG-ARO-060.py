ITEM = '/home/claude/chem5/out/items/ORG-ARO-060/item.json'
# 検証スクリプト（グループ g09）  python3 -I この.py
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
assert np_ == na, ('小問の数', np_, na)
assert [int(m.group(2)) for m in re.finditer(r'(^|\s)\((\d+)\)', it['answer'])] == list(range(1, np_+1))

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import Descriptors
cands = {'アセチルサリチル酸':'CC(=O)Oc1ccccc1C(=O)O','アセトアミノフェン':'CC(=O)Nc1ccc(O)cc1','フェナセチン':'CCOc1ccc(NC(C)=O)cc1','アセトアニリド':'CC(=O)Nc1ccccc1','スルファニルアミド':'Nc1ccc(cc1)S(N)(=O)=O'}
A = {'H':1.0,'C':12,'N':14,'O':16,'S':32}
import re as _re
hit = []
for k,s in cands.items():
    f = CalcMolFormula(Chem.MolFromSmiles(s))
    d = {e:int(n or 1) for e,n in _re.findall(r'([A-Z][a-z]?)(\d*)', f)}
    M = sum(A[e]*n for e,n in d.items())
    pc = {e: round(100*A[e]*d.get(e,0)/M,1) for e in 'CHN'}
    print(k, f, M, pc)
    if pc == {'C':67.0,'H':7.3,'N':7.8}: hit.append(k)
assert hit == ['フェナセチン']
# 官能基による実験結果の判定
pat = {'COOH':'[CX3](=O)[OX2H1]','phenolOH':'c[OX2H1]','ArNH2':'c[NX3H2]','amide':'C(=O)[NX3H1]c'}
res = {}
for k,s in cands.items():
    m = Chem.MolFromSmiles(s)
    res[k] = {p: m.HasSubstructMatch(Chem.MolFromSmarts(q)) for p,q in pat.items()}
print(res)
assert [k for k in res if res[k]['COOH']] == ['アセチルサリチル酸']
assert [k for k in res if res[k]['phenolOH']] == ['アセトアミノフェン']
assert [k for k in res if res[k]['ArNH2']] == ['スルファニルアミド']
assert sorted(k for k in res if res[k]['amide']) == sorted(['アセトアミノフェン','フェナセチン','アセトアニリド'])
for e in ['C6H5NHCOCH3 + H2O + HCl -> C6H5NH3Cl + CH3COOH','C6H4(OCOCH3)COOH + H2O -> C6H4(OH)COOH + CH3COOH']:
    assert check_equation(e) is None
print('ALL OK')
