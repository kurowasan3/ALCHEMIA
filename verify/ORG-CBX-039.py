ITEM = '/home/claude/chem6/out/items/ORG-CBX-039/item.json'
# 検証スクリプト（グループ g08）  python3 -I この.py
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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if should_check(b):
            r = check_equation(clean_ce(b))
            assert r is None, (b, r)
            n_eq += 1
            print('式OK:', b)
print('検査した式の数:', n_eq)

# 小問数と答えの数の一致
n_items = it['problem'].count('\\item')
n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))


import sympy as sp
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
n=sp.symbols('n'); M=12*n+1*(2*n-2)+32; assert sp.expand(M-(14*n+30))==0
assert sp.solve(sp.Eq(M,86),n)==[4]
assert abs(0.258/(0.100*30.0/1000)-86)<1e-9 and abs(0.430/(1.00*5.00/1000)-86)<1e-9
assert round(2.0/86*100,2)==2.33 and round(2.0/0.0233)==86
for s in ['C/C=C/C(=O)O','C=CCOC=O','C=CC(=O)OC','C=CCC(=O)O','C=C(C)C(=O)O']:
    assert CalcMolFormula(Chem.MolFromSmiles(s))=='C4H6O2'
# cis-trans を持つ酸
def has_ct(s):
    m=Chem.MolFromSmiles(s); from rdkit.Chem import FindPotentialStereo
    return any(si.type==Chem.StereoType.Bond_Double for si in Chem.FindPotentialStereo(m))
acids={'C=CCC(=O)O':False,'CC=CC(=O)O':True,'C=C(C)C(=O)O':False}
for s,v in acids.items(): print(s,has_ct(s)); assert has_ct(s)==v
# C3H5OH 鎖式アルコール（エノール除く）
cands=['C=CCO','CC=CO','C=C(C)O']
enol=Chem.MolFromSmarts('C=C[OX2H]')
ok=[s for s in cands if not Chem.MolFromSmiles(s).HasSubstructMatch(enol)]; print(ok); assert ok==['C=CCO']

print('全検証OK:', ITEM)
