# 検証スクリプト（x1 g03）  python3 -I この.py
import sys, json, re
from math import log10, floor, sqrt, pi, exp, log
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x1/out/items/ORG-CBX-068/item.json'
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
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
c=17.60*12/44; h=9.00*2/18; o=7.40-c-h
assert abs(c/12-0.4)<1e-9 and abs(h-1.0)<1e-9 and abs(o/16-0.1)<1e-9
c2=17.60*12/44; h2=7.20*2/18; o2=7.20-c2-h2
assert abs(c2/12/(o2/16)-4)<1e-9 and abs(h2/(o2/16)-8)<1e-9
assert CalcMolFormula(Chem.MolFromSmiles('CC(O)CC'))=='C4H10O'
es=['CCCC(=O)OC(C)CC','CC(C)C(=O)OC(C)CC']
assert all(CalcMolFormula(Chem.MolFromSmiles(x))=='C8H16O2' for x in es)
assert sum(len({Chem.MolToSmiles(m) for m in EnumerateStereoisomers(Chem.MolFromSmiles(x))}) for x in es)==4
assert sig(10.1/144/0.100*100,2)==70
print('ALL OK')
