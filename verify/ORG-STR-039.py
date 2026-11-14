# 検証スクリプト（グループ g03）  python3 -I この.py

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

def should_check(eq):
    if re.search(r'boxed|□|\?|？|…|\\square|\\text|\\mathrm|\$', eq):
        return False
    if re.search(r'(^|[\s(+-])R(?=[-\s(]|$)|_\{?n|(^|\s)n\s*[A-Z(]|[A-Z]\d*n\b', eq):
        return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1:
        return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None

def check_all_ce(it):
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if should_check(b):
                r = check_equation(clean_ce(b))
                assert r is None, (b, r)
                print('式OK:', b); n += 1
    print('検査した反応式:', n)

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
it = json.load(open('/home/claude/chem5/out/items/ORG-STR-039/item.json', encoding='utf-8'))
check_all_ce(it)
def nC(s): return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
prod = {'ア': 'OC(=O)CCC(=O)O', 'イ': 'OC(=O)C(O)CC(=O)O', 'ウ': 'OC(=O)C(Br)CC(=O)O', 'エ': 'OC(=O)C(Br)C(Br)C(=O)O',
        '4a': 'BrCC(C(=O)O)C(=O)O', '4b': 'CC(Br)(C(=O)O)C(=O)O'}
res = {k: nC(v) for k, v in prod.items()}; print(res)
assert res == {'ア': 0, 'イ': 1, 'ウ': 1, 'エ': 2, '4a': 0, '4b': 0}
iso = {Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(prod['エ']))}
def mirror(s): return Chem.MolToSmiles(Chem.MolFromSmiles(s.replace('@@', '!').replace('@', '@@').replace('!', '@')))
meso = [s for s in iso if mirror(s) == s]; chiral = [s for s in iso if mirror(s) != s]
print(len(iso), meso, chiral)
assert len(iso) == 3 and len(meso) == 1 and len(chiral) == 2 and mirror(chiral[0]) == chiral[1]
# 図の立体式がこの 3 種と一致すること
fig = ['O=C(O)[C@H](Br)[C@@H](Br)C(=O)O', 'O=C(O)[C@@H](Br)[C@H](Br)C(=O)O', 'O=C(O)[C@@H](Br)[C@@H](Br)C(=O)O']
assert {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in fig} == iso and mirror(fig[2]) == Chem.MolToSmiles(Chem.MolFromSmiles(fig[2]))
print('すべて一致')
