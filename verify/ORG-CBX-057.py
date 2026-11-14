# 検証スクリプト（東工大 g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/tokodai/out/items/ORG-CBX-057/item.json'
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

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans <= 1), (nprob, nans)
print('小問数', nprob, '答え', nans)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
# ---- 個別の検証 ----

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
acids = {1:['[H]C(=O)O'], 2:['CC(=O)O'], 3:['CCC(=O)O'], 4:['CCCC(=O)O','CC(C)C(=O)O']}
alcs = {1:['CO'], 2:['CCO'], 3:['CCCO','CC(C)O'], 4:['CCCCO','CC(O)CC','CC(C)CO','CC(C)(C)O']}
esters = {}
for ca in acids:
    for acid in acids[ca]:
        for alc in alcs.get(5-ca, []):
            ac = acid.replace('[H]','')[:-1]  # drop terminal O
            m = Chem.MolFromSmiles(alc)
            rxn = Chem.MolFromSmiles(ac + 'O' )
            # build ester SMILES: acyl + O-R ; alcohol written with O last or inside
            from rdkit.Chem import AllChem
            r = AllChem.ReactionFromSmarts('[C:1](=[O:2])[OH].[OH:3][C:4]>>[C:1](=[O:2])[O:3][C:4]')
            p = r.RunReactants((Chem.MolFromSmiles(acid.replace('[H]','')), m))[0][0]
            Chem.SanitizeMol(p); s = Chem.MolToSmiles(p)
            assert CalcMolFormula(p) == 'C5H10O2'
            esters[s] = alc
print(len(esters), esters); assert len(esters) == 9
iodo = Chem.MolFromSmarts('[CH3][CX4;H1,H2][OX2H]')
c1 = sum(1 for a in esters.values() if Chem.MolFromSmiles(a).HasSubstructMatch(iodo)); assert c1 == 3, c1
c2 = sum(1 for s in esters if Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True)); assert c2 == 1, c2
c3 = sum(1 for a in esters.values() if a == 'CO'); assert c3 == 2
tert = Chem.MolFromSmarts('[CX4H0]([#6])([#6])([#6])[OX2H]')
t = [s for s, a in esters.items() if Chem.MolFromSmiles(a).HasSubstructMatch(tert)]; assert len(t) == 1; print(t)

print('ALL OK')
