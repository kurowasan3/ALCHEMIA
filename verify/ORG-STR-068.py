ITEM = '/home/claude/chemx/x1/outh/items/ORG-STR-068/item.json'
# 検証スクリプト（x1 g01）  python3 -I このファイル
import sys, json, re, math
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    if x == 0: return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))

def chk(name, val, expect, n):
    v = sig(val, n)
    assert abs(v - expect) <= 1e-9 * max(1, abs(expect)), (name, val, v, expect)
    print(f'OK {name}: {val:.6g} -> {v} (期待 {expect})')

nq = sum(seg.count('\\item') for seg in re.findall(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', it['problem'], re.S))
na = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nq == na or (nq == 0 and na <= 1), ('小問数', nq, na)
print('小問数OK', nq)

# ---- 個別の検証 ----

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def mirror(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('@@','!').replace('@','@@').replace('!','@')))
SS=Chem.MolToSmiles(Chem.MolFromSmiles('CC[C@H](C)OC(=O)CCCC(=O)O[C@@H](C)CC'))
SR=Chem.MolToSmiles(Chem.MolFromSmiles('CC[C@H](C)OC(=O)CCCC(=O)O[C@H](C)CC'))
assert CalcMolFormula(Chem.MolFromSmiles(SS))=='C13H24O4'
assert mirror(SR)==SR and mirror(SS)!=SS and SS!=SR
from fractions import Fraction as F
s,r=F(3,4),F(1,4)
meso=2*s*r; chk('meso %',float(meso)*100,38,2)
ee=(s*s-r*r)/(s*s+r*r); chk('ee %',float(ee)*100,80,2)
e=F(1,2); assert ee==2*e/(1+e*e)

print('ALL OK')
