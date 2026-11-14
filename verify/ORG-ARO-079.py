# 検証スクリプト ORG-ARO-079  python3 -I この.py
ITEM = '/home/claude/chemx/kyodai/outh/items/ORG-ARO-079/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    for e in all_equations():
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK', e)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    if n_items: assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()


from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import combinations
from collections import Counter
cs = lambda m: Chem.MolToSmiles(m)
base = Chem.MolFromSmiles('c1ccccc1Oc1ccccc1')
Hpos = [a.GetIdx() for a in base.GetAtoms() if a.GetSymbol()=='C' and a.GetTotalNumHs()==1]
def add(idxs, n):
    rw = Chem.RWMol(base)
    for i in idxs:
        prev = i
        for _ in range(n):
            k = rw.AddAtom(Chem.Atom(6)); rw.AddBond(prev, k, Chem.BondType.SINGLE); prev = k
    m = rw.GetMol(); Chem.SanitizeMol(m); return cs(m)
iso = {add([i],2) for i in Hpos} | {add([i,j],1) for i,j in combinations(Hpos,2)}
assert len(iso) == 15 and all(CalcMolFormula(Chem.MolFromSmiles(s))=='C14H14O' for s in iso)
def products(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    O = [a for a in m.GetAtoms() if a.GetSymbol()=='O'][0]
    for nb in O.GetNeighbors():
        rw = Chem.RWMol(m); rw.RemoveBond(O.GetIdx(), nb.GetIdx()); mm = rw.GetMol(); Chem.SanitizeMol(mm)
        out |= {cs(f) for f in Chem.GetMolFrags(mm, asMols=True)}
    return out
cnt = Counter(len(products(s)) for s in iso); print(cnt)
assert cnt == {2:3, 3:3, 4:9}
def phenols(ps): return [p for p in ps if Chem.MolFromSmiles(p).HasSubstructMatch(Chem.MolFromSmarts('c[OX2H]'))]
def brominate(smi):
    m = Chem.MolFromSmiles(smi); rw = Chem.RWMol(m)
    ri = m.GetRingInfo()
    cOH = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=='C' and any(n.GetSymbol()=='O' for n in a.GetNeighbors())][0]
    dm = Chem.GetDistanceMatrix(m)
    ring = [r for r in ri.AtomRings() if cOH in r][0]
    targets = [i for i in ring if dm[cOH][i] in (1,3) and m.GetAtomWithIdx(i).GetTotalNumHs()==1]
    for i in targets:
        k = rw.AddAtom(Chem.Atom(35)); rw.AddBond(i, k, Chem.BondType.SINGLE)
    mm = rw.GetMol(); Chem.SanitizeMol(mm); return CalcMolFormula(mm)
def nCl(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    for a in m.GetAtoms():
        if a.GetIsAromatic() and a.GetTotalNumHs()==1:
            rw = Chem.RWMol(m); k = rw.AddAtom(Chem.Atom(17)); rw.AddBond(a.GetIdx(), k, Chem.BondType.SINGLE)
            mm = rw.GetMol(); Chem.SanitizeMol(mm); out.add(cs(mm))
    return len(out)
# (3)
P = [s for s in iso if len(products(s))==2 and all(brominate(p)=='C7H5Br3O' for p in phenols(products(s)))]
print('P', P); assert P == [cs(Chem.MolFromSmiles('Cc1cccc(Oc2cccc(C)c2)c1'))]
# (4)
Q = []
for s in iso:
    ps = products(s)
    if len(ps)!=4: continue
    ph = phenols(ps); T = [p for p in ph if CalcMolFormula(Chem.MolFromSmiles(p))=='C8H10O']
    if cs(Chem.MolFromSmiles('Oc1ccccc1')) in ph and T and brominate(T[0])=='C8H9BrO' and nCl(T[0])==2: Q.append(s)
print('Q', Q); assert Q == [cs(Chem.MolFromSmiles('Cc1cccc(C)c1Oc1ccccc1'))]
# (5)
U = []
for s in iso:
    ps = products(s)
    if len(ps)!=3: continue
    ph = phenols(ps); f = {p: brominate(p) for p in ph}
    if sorted(f.values()) == ['C7H5Br3O','C7H6Br2O'] and nCl([p for p in ph if f[p]=='C7H6Br2O'][0])==2: U.append(s)
print('U', U); assert U == [cs(Chem.MolFromSmiles('Cc1ccc(Oc2cccc(C)c2)cc1'))]
print('OK')
