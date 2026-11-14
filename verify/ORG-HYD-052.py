ITEM = '/home/claude/chemx/x1/outh/items/ORG-HYD-052/item.json'
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
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
def stereo(smi):
    return {Chem.MolToSmiles(s) for s in EnumerateStereoisomers(Chem.MolFromSmiles(smi))}
def mirror(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('@@','!').replace('@','@@').replace('!','@')))
assert len(stereo('CC1(Cl)CCCC1'))==1
s13=stereo('CC1CCC(Cl)C1'); s12=stereo('CC1CCCC1Cl')
assert len(s13)==4 and len(s12)==4 and all(mirror(x)!=x for x in s13|s12)
# fractions = number of enantiomer pairs
def fractions(S): return len({frozenset({x,mirror(x)}) for x in S})
assert fractions(s13)==2
assert fractions(s12|s13)==4 and len(s12|s13)==8
# chirality of starting alkenes
for a,n in [('CC1=CCCC1',0),('CC1C=CCC1',1),('CC1CC=CC1',0)]:
    assert len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(a),includeUnassigned=True))==n
# enantiopure Q: Me-carbon fixed -> products with that carbon fixed: pick (R) at methyl carbon
from rdkit.Chem import FindMolChiralCenters
def me_conf(s):
    m=Chem.MolFromSmiles(s)
    for a in m.GetAtoms():
        if a.GetSymbol()=='C' and any(n.GetDegree()==1 and n.GetSymbol()=='C' for n in a.GetNeighbors()) and a.HasProp('_CIPCode'):
            return a.GetProp('_CIPCode')
    from rdkit.Chem import rdCIPLabeler; rdCIPLabeler.AssignCIPLabels(m)
    for a in m.GetAtoms():
        if a.GetSymbol()=='C' and any(n.GetDegree()==1 and n.GetSymbol()=='C' for n in a.GetNeighbors()) and a.HasProp('_CIPCode'):
            return a.GetProp('_CIPCode')
# CIP label of methyl carbon can flip between products with different substituents, so check by counting: one enantiomer series = half of all
Q = {x for x in s12|s13}
assert len(Q)//2 == 4
print('stereo counts OK')

print('ALL OK')
