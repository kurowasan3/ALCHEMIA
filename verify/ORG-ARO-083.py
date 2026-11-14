# 検証スクリプト（東工大 g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/out/items/ORG-ARO-083/item.json'
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
from rdkit.Chem import Descriptors, rdMolDescriptors
from itertools import combinations
def can(m): return Chem.MolToSmiles(m)
def mono_cl(smi):
    m = Chem.MolFromSmiles(smi); out = set()
    for a in m.GetAtoms():
        if a.GetIsAromatic() and a.GetTotalNumHs() == 1:
            rw = Chem.RWMol(m); i = rw.AddAtom(Chem.Atom(17)); rw.AddBond(a.GetIdx(), i, Chem.BondType.SINGLE)
            mm = rw.GetMol(); Chem.SanitizeMol(mm); out.add(can(mm))
    return len(out)
# enumerate dihydroxybenzoic acids: ring positions 0..5, COOH at 0
acids = set()
for pos in combinations(range(1, 6), 2):
    subs = ['C(=O)O' if i == 0 else ('O' if i in pos else '') for i in range(6)]
    s = 'c1(%s)' % subs[0]
    ring = ''.join('c(%s)' % subs[i] if subs[i] else 'c' for i in range(1, 6))
    acids.add(can(Chem.MolFromSmiles(s + ring + '1')))
print('acids', len(acids)); assert len(acids) == 6
prod = {}
for a in acids:
    m = Chem.MolFromSmiles(a)
    M = Descriptors.MolWt(m); assert abs(M - 154.12) < 0.05
    assert round(44.01 / M * 100, 1) == 28.6, 44.01/M
    patt = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    q = can(Chem.RemoveHs(Chem.DeleteSubstructs(m, patt)))
    q = can(Chem.MolFromSmiles(q))
    prod.setdefault(q, []).append(a)
names = {can(Chem.MolFromSmiles('Oc1ccccc1O')): 'cat', can(Chem.MolFromSmiles('Oc1cccc(O)c1')): 'res', can(Chem.MolFromSmiles('Oc1ccc(O)cc1')): 'hyd'}
cnt = {names[k]: len(v) for k, v in prod.items()}
print(cnt); assert cnt == {'cat': 2, 'res': 3, 'hyd': 1}
assert mono_cl('Oc1ccc(O)cc1') == 1 and mono_cl('Oc1ccccc1O') == 2 and mono_cl('Oc1cccc(O)c1') == 3
assert mono_cl('OC(=O)c1cc(O)ccc1O') == 3
print('Q MW', Descriptors.MolWt(Chem.MolFromSmiles('Oc1ccc(O)cc1')), rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles('Oc1ccc(O)cc1')))
print('OK')
