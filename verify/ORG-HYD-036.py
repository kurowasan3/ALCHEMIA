ITEM = '/home/claude/chem6/out/items/ORG-HYD-036/item.json'
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


from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
can=lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
add=AllChem.ReactionFromSmarts('[C:1]#[C:2]>>[C:1]([H])=[C:2]Cl')
def mono(s):
    out=set()
    for ps in add.RunReactants((Chem.MolFromSmiles(s),)):
        p=Chem.RemoveHs(ps[0]); Chem.SanitizeMol(p)
        for iso in EnumerateStereoisomers(p): out.add(Chem.MolToSmiles(iso))
    return out
m1=mono('C#CCC'); m2=mono('CC#CC'); print(m1,m2)
assert len(m1)==3 and len(m2)==2
gem=can('CC(Cl)(Cl)CC')
assert can('C=C(Cl)CC'.replace('C=C(Cl)','CC(Cl)(Cl)'))==gem
assert can('CC(Cl)=CC'.replace('C(Cl)=CC','C(Cl)(Cl)CC'))==gem
for s in [gem, can('ClC=C(Cl)CC')]:
    m=Chem.MolFromSmiles(s)
    n=len(Chem.FindMolChiralCenters(m,includeUnassigned=True,useLegacyImplementation=False)); print(s,n); assert n==0
for e in ['C4H6 + HCl -> C4H7Cl','C4H7Cl + HCl -> C4H8Cl2','C4H6 + Cl2 -> C4H6Cl2']:
    assert check_equation(e) is None

print('全検証OK:', ITEM)
