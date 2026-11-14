ITEM = '/home/claude/chem2/out/items/ORG-STR-018/item.json'
# 検証スクリプト（グループ g08）  python3 -I この.py
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

def check_all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return eqs

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def check_parts():
    n_items = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, n_items + 1)], (n_items, nums)
    print('小問の数OK:', n_items)


check_parts()
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import itertools
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
mirror = lambda s: can(s.replace('@@', '!').replace('@', '@@').replace('!', '@'))
CELLS = [('C[C@@H](O)[C@@H](C)Cl', 'ア'), ('C[C@H](O)[C@H](C)Cl', 'イ'), ('C[C@@H](O)[C@H](C)Cl', 'ウ'), ('Cl[C@H](C)[C@H](O)C', 'エ'), ('Cl[C@H](C)[C@@H](O)C', 'オ')]

S = {lab: can(smi) for smi, lab in CELLS}
same = sorted(a + b for a, b in itertools.combinations(S, 2) if S[a] == S[b])
enan = sorted(a + b for a, b in itertools.combinations(S, 2) if mirror(S[a]) == S[b] and S[a] != S[b])
print('同一', same, '鏡像', enan); assert same == ['アエ'] and enan == ['アイ', 'イエ', 'ウオ']
# 図の SVG が CELLS から作られたことを確認
svg = open(ITEM.replace('item.json', 'problem.svg'), encoding='utf-8').read()
for _, lab in CELLS: assert '>' + lab + '<' in svg
def st(s): return sorted({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(Chem.MolFromSmiles(s), options=StereoEnumerationOptions(unique=True))})
assert len(st('CC(O)C(C)Cl')) == 4
bd = st('CC(O)C(C)O'); assert len(bd) == 3 and sum(mirror(x) == x for x in bd) == 1

print('全検証OK:', ITEM)
