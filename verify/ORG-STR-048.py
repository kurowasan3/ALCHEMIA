ITEM = '/home/claude/chem6/out/items/ORG-STR-048/item.json'
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
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
import itertools
pairs=list(itertools.combinations(range(5),2))
found=set()
for orders in itertools.product((0,1,2),repeat=10):
    deg=[0]*5; adj={i:set() for i in range(5)}
    for (i,j),o in zip(pairs,orders):
        if o: deg[i]+=o; deg[j]+=o; adj[i].add(j); adj[j].add(i)
    if max(deg)>4 or sum(4-d for d in deg)!=10: continue
    seen={0}; st=[0]
    while st:
        u=st.pop()
        for v in adj[u]:
            if v not in seen: seen.add(v); st.append(v)
    if len(seen)<5: continue
    rw=Chem.RWMol()
    for _ in range(5): rw.AddAtom(Chem.Atom(6))
    for (i,j),o in zip(pairs,orders):
        if o: rw.AddBond(i,j,Chem.BondType.SINGLE if o==1 else Chem.BondType.DOUBLE)
    try:
        m=rw.GetMol(); Chem.SanitizeMol(m); found.add(Chem.MolToSmiles(m))
    except Exception: pass
print(len(found), sorted(found)); assert len(found)==10
ring=[s for s in found if Chem.MolFromSmiles(s).GetRingInfo().NumRings()>0]; assert len(ring)==5
allst=set()
for s in found:
    for iso in EnumerateStereoisomers(Chem.MolFromSmiles(s)):
        allst.add(Chem.MolToSmiles(iso))
print(len(allst), sorted(allst)); assert len(allst)==13
def mirror(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('@@','!').replace('@','@@').replace('!','@')))
cis='C[C@@H]1C[C@@H]1C'; assert Chem.MolToSmiles(Chem.MolFromSmiles(cis))==mirror(cis)
trans='C[C@@H]1C[C@H]1C'; assert Chem.MolToSmiles(Chem.MolFromSmiles(trans))!=mirror(trans)
# (2) 水素化で2-メチルブタン
alk=[s for s in found if s not in ring]
def hyd(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('=','')))
ipent=Chem.MolToSmiles(Chem.MolFromSmiles('CC(C)CC'))
r2=sorted(s for s in alk if hyd(s)==ipent); print(r2); assert len(r2)==3
# (3) 臭素付加物の不斉炭素数
from rdkit.Chem import AllChem
rxn=AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1](Br)[C:2]Br')
nc={}
for s in alk:
    p=rxn.RunReactants((Chem.MolFromSmiles(s),))[0][0]; Chem.SanitizeMol(p)
    nc[s]=len(Chem.FindMolChiralCenters(p,includeUnassigned=True,useLegacyImplementation=False))
print(nc); assert [s for s,v in nc.items() if v==2]==[Chem.MolToSmiles(Chem.MolFromSmiles('CC=CCC'))]

print('全検証OK:', ITEM)
