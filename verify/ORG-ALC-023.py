# 検証スクリプト（グループ g02）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-ALC-023/item.json'
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
        if ' + ' not in b and not re.search(r'(^|\s)\d+\s*[A-Z(\[]', b): continue
        r = check_equation(b); assert r is None, (b, r); n_eq += 1; print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
from math import log10, floor
def sig(x, n): return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import product
# C5H12O の鎖式アルコールを列挙（炭素骨格 + OH の位置）
skel = ['CCCCC', 'CC(C)CC', 'CC(C)(C)C']
alcs = set()
for s in skel:
    m = Chem.MolFromSmiles(s)
    for a in m.GetAtoms():
        if a.GetTotalNumHs() > 0:
            em = Chem.RWMol(m); i = em.AddAtom(Chem.Atom(8)); em.AddBond(a.GetIdx(), i, Chem.BondType.SINGLE)
            alcs.add(Chem.MolToSmiles(em))
assert len(alcs) == 8, alcs
def cls(smi):
    m = Chem.MolFromSmiles(smi); o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    c = o.GetNeighbors()[0]; return sum(1 for n in c.GetNeighbors() if n.GetSymbol() == 'C')
sec = sorted(a for a in alcs if cls(a) == 2)
assert len(sec) == 3, sec
def ketone(smi):
    m = Chem.RWMol(Chem.MolFromSmiles(smi)); o = [a for a in m.GetAtoms() if a.GetSymbol() == 'O'][0]
    b = m.GetBondBetweenAtoms(o.GetIdx(), o.GetNeighbors()[0].GetIdx()); b.SetBondType(Chem.BondType.DOUBLE)
    Chem.SanitizeMol(m); return Chem.MolToSmiles(m)
noiod = [s for s in sec if not Chem.MolFromSmiles(ketone(s)).HasSubstructMatch(Chem.MolFromSmarts('[CH3]C(=O)'))]
assert [Chem.MolToSmiles(Chem.MolFromSmiles(x)) for x in noiod] == [Chem.MolToSmiles(Chem.MolFromSmiles('CCC(O)CC'))], noiod
# C4H10O の級数
c4 = {'1-ブタノール':'CCCCO','2-ブタノール':'CCC(C)O','2-メチル-1-プロパノール':'CC(C)CO','2-メチル-2-プロパノール':'CC(C)(C)O'}
print({k: cls(v) for k, v in c4.items()})
assert [cls(v) for v in c4.values()] == [1, 2, 1, 3]
print('P = 3-ペンタノール  OK')
