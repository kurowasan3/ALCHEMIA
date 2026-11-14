# 検証スクリプト（東工大 g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/outh/items/ORG-STR-063/item.json'
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

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans <= 1), (nprob, nans)
print('小問数', nprob, '答え', nans)
for k in ['src_no', 'src_label', 'field', 'theme', 'difficulty', 'time_min', 'format', 'problem', 'answer', 'solution', 'core', 'tags', 'verified', 'verify_note', 'problem_fig', 'solution_fig']:
    assert k in it, k
# ---- 個別の検証 ----

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
m = Chem.AddHs(Chem.MolFromSmiles('CC(=O)Nc1ccccc1'))
E = sum(a.GetAtomicNum() for a in m.GetAtoms())
val = {6: 4, 7: 3, 8: 2, 1: 1}
S = sum(val[a.GetAtomicNum()] for a in m.GetAtoms())
Chem.Kekulize(m, clearAromaticFlags=True)
B = sum(b.GetBondTypeAsDouble() for b in m.GetBonds())
print(E, S, B); assert (E, S, B) == (72, 46, 23)
sols = [(x, y, z) for x in range(1, 30) for y in range(1, 70) for z in range(0, 10) if 6*x + y + 8*z == 72 and 4*x + y + 2*z == 46 and y <= 2*x + 2]
print(sols); assert sorted(sols) == [(7, 14, 2), (10, 4, 1)]
buts = {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in ['CCCCO', 'CCC(C)O', 'CC(C)CO', 'CC(C)(C)O']}
assert len(buts) == 4 and all(rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(s)) == 'C4H10O' for s in buts)
assert round(Descriptors.MolWt(Chem.MolFromSmiles('CCC(=O)O'))) == 74 and round(Descriptors.MolWt(Chem.MolFromSmiles('CCCCO'))) == 74
# equal MW (integer masses): acid CnH2n+1COOH = 14n+46, alcohol CmH2m+1OH = 14m+18, n+m=6
pairs = [(n, 6 - n) for n in range(0, 7) if 14*n + 46 == 14*(6 - n) + 18]
print(pairs); assert pairs == [(2, 4)]
C = Chem.MolFromSmiles('CCC(=O)OC(C)CC')
assert rdMolDescriptors.CalcMolFormula(C) == 'C7H14O2'
cc = Chem.FindMolChiralCenters(C, includeUnassigned=True, useLegacyImplementation=False); print(cc); assert len(cc) == 1
# C's own electrons and valence sum
Ch = Chem.AddHs(C)
assert sum(a.GetAtomicNum() for a in Ch.GetAtoms()) == 72 and sum(val[a.GetAtomicNum()] for a in Ch.GetAtoms()) == 46
print('OK')
