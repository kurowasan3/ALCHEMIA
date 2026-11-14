ITEM = '/home/claude/chemx/x1/out/items/POL-PRT-070/item.json'
# 検証スクリプト（x1 g04）  python3 -I この.py
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
assert nq == na or (nq == 0 and na == 0), ('小問数', nq, na)
print('小問数OK', nq)

# ---- 個別の検証 ----

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem import FindMolChiralCenters
for e in ['CH3CHO + NH3 + HCN -> CH3CH(NH2)CN + H2O','CH3CH(NH2)CN + 2H2O + HCl -> CH3CH(NH2)COOH + NH4Cl']:
    assert check_equation(e) is None, e
def strecker(ald_smiles):
    m = Chem.MolFromSmiles(ald_smiles)
    rxn = __import__('rdkit.Chem.AllChem', fromlist=['x']).ReactionFromSmarts('[CX3:1]=[O]>>[C:1](N)C(=O)O')
    p = rxn.RunReactants((m,))[0][0]; Chem.SanitizeMol(p); return Chem.CanonSmiles(Chem.MolToSmiles(p))
cases = {'CC(C)C=O':'CC(C)C(N)C(=O)O','O=CCc1ccccc1':'NC(Cc1ccccc1)C(=O)O','C=O':'NCC(=O)O','CC(C)=O':'CC(C)(N)C(=O)O'}
nch = []
for a,t in cases.items():
    p = strecker(a); assert p == Chem.CanonSmiles(t), (a,p)
    nch.append(len(FindMolChiralCenters(Chem.MolFromSmiles(p), includeUnassigned=True, useLegacyImplementation=False)))
assert nch == [1,1,0,0], nch
print('Strecker 生成物OK', nch)
vf = Chem.MolFromSmiles('CC(C)C(N)C(=O)NC(Cc1ccccc1)C(=O)O')
isos = [Chem.MolToSmiles(x) for x in EnumerateStereoisomers(vf)]
assert len(set(isos)) == 4
def mirror(s): return Chem.CanonSmiles(s.replace('@@','!').replace('@','@@').replace('!','@'))
pairs = {frozenset((Chem.CanonSmiles(s), mirror(s))) for s in isos}
assert len(pairs) == 2 and all(len(p)==2 for p in pairs)
print('Val-Phe 立体異性体 4, 鏡像の組 2')

print('ALL OK')
