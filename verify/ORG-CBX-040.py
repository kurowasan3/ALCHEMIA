ITEM = '/home/claude/chem6/out/items/ORG-CBX-040/item.json'
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


C=9.63*12/44; H=3.38*2.0/18; Ox=5.00-C-H
r=[C/12,H/1.0,Ox/16]; mn=min(r); rr=[x/mn for x in r]; print(rr)
assert [round(x*4) for x in rr]==[7,12,4] and all(abs(x*4-round(x*4))<0.1 for x in rr)
assert 7*12+12+64==160
nNaOH=0.100*50.0/1000; M_per_x=0.400/nNaOH; print(M_per_x); assert abs(M_per_x-80)<1e-9
for e in ['(CH2)5(COO)2Ca -> (CH2)5CO + CaCO3','CaC7H10O4 -> C6H10O + CaCO3']:
    assert check_equation(e) is None
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
rxn=AllChem.ReactionFromSmarts('([CX3:1](=[O:5])[OX2H1].[C:3][CX3:4](=O)[OX2H1])>>[C:1](=[O:5])[C:3]')
def prods(s):
    ps=set()
    for p in rxn.RunReactants((Chem.MolFromSmiles(s),)):
        m=p[0]; Chem.SanitizeMol(m); ps.add(Chem.MolToSmiles(m))
    return ps
can=lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
assert CalcMolFormula(Chem.MolFromSmiles('OC(=O)CCCCCC(=O)O'))=='C7H12O4'
assert prods('OC(=O)CCCCCC(=O)O')=={can('O=C1CCCCC1')}
assert prods('OC(=O)CCCCC(=O)O')=={can('O=C1CCCC1')}
assert prods('OC(=O)C(C)CCCC(=O)O')=={can('CC1CCCC1=O')}
assert prods('OC(=O)CC(C)CCC(=O)O')=={can('CC1CCC(=O)C1')}
assert CalcMolFormula(Chem.MolFromSmiles('O=C1CCCCC1'))=='C6H10O'
for m in range(3,7):
    s='OC(=O)'+'C'*m+'C(=O)O'; p=list(prods(s))[0]; pm=Chem.MolFromSmiles(p); ri=pm.GetRingInfo()
    assert ri.AtomRings() and len(ri.AtomRings()[0])==m+1
print('ring size m+1 OK')

print('全検証OK:', ITEM)
