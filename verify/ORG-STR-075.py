# 検証スクリプト（グループ g05）  python3 -I このファイル
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemk/out/items/ORG-STR-075/item.json'
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
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from collections import Counter
iso = {'ヘプタン':'CCCCCCC','2-メチルヘキサン':'CC(C)CCCC','3-メチルヘキサン':'CCC(C)CCC','3-エチルペンタン':'CCC(CC)CC',
 '2,2-ジメチルペンタン':'CC(C)(C)CCC','3,3-ジメチルペンタン':'CCC(C)(C)CC','2,3-ジメチルペンタン':'CC(C)C(C)CC',
 '2,4-ジメチルペンタン':'CC(C)CC(C)C','2,2,3-トリメチルブタン':'CC(C)(C)C(C)C'}
can = set()
for n,s in iso.items():
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m)=='C7H16'; can.add(Chem.MolToSmiles(m))
assert len(can)==9
# 全列挙：炭素鎖に C を1個ずつ付けていく
cur={'C'}
for _ in range(6):
    nxt=set()
    for s in cur:
        m=Chem.MolFromSmiles(s)
        for a in m.GetAtoms():
            if a.GetDegree()<4:
                rw=Chem.RWMol(m); i=rw.AddAtom(Chem.Atom(6)); rw.AddBond(a.GetIdx(),i,Chem.BondType.SINGLE)
                nxt.add(Chem.MolToSmiles(rw))
    cur=nxt
assert cur==can, (len(cur))
types=Counter()
for n,s in iso.items():
    m=Chem.MolFromSmiles(s)
    deg=[a.GetDegree() for a in m.GetAtoms()]
    t=deg.count(3); q=deg.count(4); me=deg.count(1)
    assert me==2+t+2*q
    types[(t,q)]+=1
print(types)
assert types=={(0,0):1,(1,0):3,(0,1):2,(2,0):2,(1,1):1}
def nmono(s):
    m=Chem.AddHs(Chem.MolFromSmiles(s)); out=set()
    for a in m.GetAtoms():
        if a.GetSymbol()=='H':
            rw=Chem.RWMol(m); rw.GetAtomWithIdx(a.GetIdx()).SetAtomicNum(17)
            out.add(Chem.MolToSmiles(Chem.RemoveHs(rw.GetMol())))
    return len(out)
five=[n for n,s in iso.items() if [a.GetDegree() for a in Chem.MolFromSmiles(s).GetAtoms()].count(1)==5]
assert five==['2,2,3-トリメチルブタン']
assert nmono(iso['2,2,3-トリメチルブタン'])==3
four_noCH=[n for n,s in iso.items() if (lambda d:(d.count(1)==4 and d.count(3)==0))([a.GetDegree() for a in Chem.MolFromSmiles(s).GetAtoms()])]
assert sorted(four_noCH)==sorted(['2,2-ジメチルペンタン','3,3-ジメチルペンタン'])
assert nmono(iso['2,2-ジメチルペンタン'])==4 and nmono(iso['3,3-ジメチルペンタン'])==3
print('異性体・分類・置換体 OK')

print('ALL OK')
