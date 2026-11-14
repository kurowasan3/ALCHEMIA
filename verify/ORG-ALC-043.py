ITEM = '/home/claude/chemk/out/items/ORG-ALC-043/item.json'
# 検証スクリプト（グループ g04）  python3 -I この.py
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
F = lambda s: CalcMolFormula(Chem.MolFromSmiles(s))
def add(*fs):
    c = Counter()
    for f in fs:
        for e, n in re.findall(r'([A-Z][a-z]?)(\d*)', f): c[e] += int(n or 1)
    return c
A = 'CCC(O)C(C)C=O'; B = 'CCC=C(C)C=O'; C = 'CC(=O)C=Cc1ccccc1'; D = 'O=C(C=Cc1ccccc1)C=Cc1ccccc1'
assert F(A) == 'C6H12O2' and F(B) == 'C6H10O' and F(C) == 'C10H10O' and F(D) == 'C17H14O'
assert len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(A), includeUnassigned=True)) == 2
assert add(F('CCC=O'), F('CCC=O')) == add(F(A))
assert add(F('CC(C)=O'), F('O=Cc1ccccc1')) == add(F(C), 'H2O')
assert add(F('CC(C)=O'), F('O=Cc1ccccc1'), F('O=Cc1ccccc1')) == add(F(D), 'H2O', 'H2O')
def alphaH(s):
    m = Chem.AddHs(Chem.MolFromSmiles(s)); n = 0
    for a in m.GetAtoms():
        if a.GetSymbol() == 'C' and any(nb.GetSymbol() == 'O' and m.GetBondBetweenAtoms(a.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 2 for nb in a.GetNeighbors()):
            for nb in a.GetNeighbors():
                if nb.GetSymbol() == 'C':
                    n += sum(1 for x in nb.GetNeighbors() if x.GetSymbol() == 'H')
    return n
cands = {'ア':'C=O','イ':'CC(C)=O','ウ':'O=Cc1ccccc1','エ':'CC(C)(C)C=O','オ':'CCC=O'}
noH = [k for k, s in cands.items() if alphaH(s) == 0]; print(noH); assert noH == ['ア','ウ','エ']

print('全検証OK:', ITEM)
