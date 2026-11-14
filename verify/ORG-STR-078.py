# 検証スクリプト（グループ g02）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation

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

def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq

def check_all_equations(it, extra=()):
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|…', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]', b): continue
            if not (' + ' in b or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', b)): continue
            e = clean_ce(b)
            r = check_equation(e)
            assert r is None, (b, r)
            print('式OK:', b); n += 1
    for e in extra:
        r = check_equation(e); assert r is None, (e, r); print('式OK(追加):', e); n += 1
    return n

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def check_parts(it):
    k = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, k + 1)], (k, nums)
    print('小問数OK:', k)
ITEM = '/home/claude/chemk/outh/items/ORG-STR-078/item.json'
it = json.load(open(ITEM, encoding='utf-8'))
check_parts(it); check_all_equations(it)
ee0 = 12.15/13.5; xs = (1 + ee0)/2; print('ee', ee0, 'xS', xs); assert round(xs*100) == 95
ee1 = -11.34/13.5; one_f = -ee1/ee0; f = 1 - one_f; xr = (1 - ee1)/2
print('f', f, 'xR', xr); assert sig(f*100, 2) == 6.7 and round(xr*100) == 92
rot2 = 13.5*ee0*one_f**2; print('2回', rot2); assert sig(rot2, 3) == 10.6
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
def side(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=7); AllChem.MMFFOptimizeMolecule(m)
    conf = m.GetConformer(); ring = [r for r in m.GetRingInfo().AtomRings() if all(m.GetAtomWithIdx(i).GetSymbol() == 'C' and not m.GetAtomWithIdx(i).GetIsAromatic() for i in r)][0]
    P = np.array([list(conf.GetAtomPosition(i)) for i in ring]); c = P.mean(0)
    n = np.linalg.svd(P - c)[2][2]
    subs = []
    for i in ring:
        for nb in m.GetAtomWithIdx(i).GetNeighbors():
            if nb.GetIdx() not in ring and nb.GetSymbol() != 'H':
                subs.append(np.dot(np.array(list(conf.GetAtomPosition(nb.GetIdx()))) - np.array(list(conf.GetAtomPosition(i))), n))
    assert len(subs) == 2
    return 'シス' if subs[0]*subs[1] > 0 else 'トランス'
cis_alc = 'C[C@H]1CC[C@@H](O)CC1'; mits = 'C[C@H]1CC[C@H](OC(=O)c2ccccc2)CC1'; fisch = 'C[C@H]1CC[C@@H](OC(=O)c2ccccc2)CC1'
r = {s: side(s) for s in (cis_alc, mits, fisch)}; print(r)
assert r[cis_alc] == 'シス' and r[mits] == 'トランス' and r[fisch] == 'シス'
for s in (cis_alc, mits):
    m = Chem.MolFromSmiles(s)
    mirror = Chem.MolFromSmiles(s.replace('@@', '!').replace('@', '@@').replace('!', '@'))
    assert Chem.MolToSmiles(m) == Chem.MolToSmiles(mirror)
    from rdkit.Chem import FindMolChiralCenters
    print(s, 'キラル中心(CIP)', FindMolChiralCenters(m, useLegacyImplementation=False))
print('鏡像と重なる（光学不活性）OK')
print('全検証OK:', ITEM)
