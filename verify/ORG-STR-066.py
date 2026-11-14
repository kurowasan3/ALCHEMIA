# 検証スクリプト（x1 g03）  python3 -I この.py
import sys, json, re
from math import log10, floor, sqrt, pi, exp, log
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x1/out/items/ORG-STR-066/item.json'
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans <= 1), (nprob, nans)
print('小問数', nprob, '答え', nans)

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
def mirror(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi.replace('@@','!').replace('@','@@').replace('!','@')))
def info(s):
    m=Chem.MolFromSmiles(s); iso=sorted({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m)})
    nc=len(Chem.FindMolChiralCenters(m,includeUnassigned=True,useLegacyImplementation=False))
    return nc,iso,[mirror(x)==x for x in iso]
a=info('ClC1(Cl)CC1'); b=info('ClC1C(Cl)C1'); c=info('BrC1C(Cl)C1')
assert a[0]==0 and len(a[1])==1
assert b[0]==2 and len(b[1])==3 and sum(b[2])==1
assert c[0]==2 and len(c[1])==4 and sum(c[2])==0
assert len(a[1])+len(b[1])==4
print('ALL OK')
