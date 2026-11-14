# 検証スクリプト（東大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/todai/outh/items/ORG-ARO-071/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

# 小問の数と答えの数の一致
nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans == 0), (nprob, nans)
print('小問数', nprob, '答え', nans)

# ---- 個別の検証 ----
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from itertools import product
def ring_subst(subs):
    # 6 員環の各位置に置換基（'' は H）を置いた全パターンを正規化して集める
    out = set()
    for pos in product(range(6), repeat=len(subs)):
        if len(set(pos)) < len(subs): continue
        at = ['c'] * 6
        for p, s in zip(pos, subs): at[p] = 'c(' + s + ')'
        smi = at[0].replace('c', 'c1', 1) + ''.join(at[1:5]) + at[5].replace('c', 'c1', 1) if False else None
        # 環を明示的に組む
        parts = []
        for i, a in enumerate(at):
            base = 'c1' if i == 0 else 'c'
            parts.append(base + (a[1:] if a != 'c' else ''))
        smi = ''.join(parts[:5]) + parts[5].replace('c', 'c1', 1)
        m = Chem.MolFromSmiles(smi)
        out.add(Chem.MolToSmiles(m))
    return out
acids = set()
for subs in (['CCC(=O)O'], ['C(C)C(=O)O'], ['CC(=O)O', 'C'], ['C(=O)O', 'CC'], ['C(=O)O', 'C', 'C']):
    acids |= ring_subst(subs)
acids = {a for a in acids if CalcMolFormula(Chem.MolFromSmiles(a)) == 'C9H10O2'}
print(len(acids)); assert len(acids) == 14
def chiral(s):
    return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False)) > 0
mono_acid = [a for a in acids if len(Chem.MolFromSmiles(a).GetSubstructMatches(Chem.MolFromSmarts('c[!c;!#1]'))) == 1]
print(mono_acid); assert len(mono_acid) == 2
ch = [a for a in acids if chiral(a)]; print(ch); assert ch == [Chem.MolToSmiles(Chem.MolFromSmiles('CC(C(=O)O)c1ccccc1'))]
alcs = set()
for subs in (['CCO'], ['C(C)O'], ['CO', 'C']):
    alcs |= ring_subst(subs)
alcs = {a for a in alcs if CalcMolFormula(Chem.MolFromSmiles(a)) == 'C8H10O'}
print(len(alcs)); assert len(alcs) == 5
cha = [a for a in alcs if chiral(a)]; assert cha == [Chem.MolToSmiles(Chem.MolFromSmiles('CC(O)c1ccccc1'))]
A = Chem.MolFromSmiles('CC(C(=O)OC(C)c1ccccc1)c1ccccc1')
assert CalcMolFormula(A) == 'C17H18O2'
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
iso = {Chem.MolToSmiles(m) for m in EnumerateStereoisomers(A)}; print(len(iso)); assert len(iso) == 4
# 元素分析・分子量
assert abs(108/150*100 - 72.0) < 0.05 and abs(10/150*100 - 6.7) < 0.05
assert abs(96/122*100 - 78.7) < 0.05 and abs(10/122*100 - 8.2) < 0.05
M = 0.150/(0.100*10.0e-3); assert abs(M - 150) < 1e-6
Mapp = 5.12*0.300/(0.256*0.0200); print(Mapp); assert abs(Mapp - 300) < 1

print('ALL OK')
