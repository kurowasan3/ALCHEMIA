# 検証スクリプト POL-SUG-029  python3 -I この.py
ITEM = '/home/claude/chem4/out/items/POL-SUG-029/item.json'
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
    """矢印が1つで完全な式（\\boxed・\\text・\\mathrm を含むものは除く）"""
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and not re.search(r'boxed|\\text|\\mathrm', b):
                eqs.add(b)
    return sorted(eqs)
def check_all_equations():
    eqs = all_equations()
    for e in eqs:
        r = check_equation(re.sub(r'(->|<=>|<-)\[[^\]]*\]', r'\1', e))
        assert r is None, (e, r)
    print('反応式', len(eqs), '本 OK:', eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem'].split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in it['problem'] else 0
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i+1) for i in range(len(nums))], nums
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

# ---- 個別の検証 ----
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
def nc(s):
    return len(Chem.FindMolChiralCenters(Chem.MolFromSmiles(s), includeUnassigned=True, useLegacyImplementation=False))
assert nc('OCC(O)C(O)C(O)C(O)C=O') == 4
assert len(list(EnumerateStereoisomers(Chem.MolFromSmiles('OCC(O)C(O)C(O)C(O)C=O')))) == 16
assert nc('OCC1OC(O)C(O)C(O)C1O') == 5
assert nc('OCC(=O)C(O)C(O)C(O)CO') == 3
assert nc('OCC1(O)OC(CO)C(O)C1O') == 4
assert nc('OCC(O)C(O)C(O)C(O)CO') == 4
can = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
# PubChem 由来の D-グルコース, D-マンノース, D-フルクトース（鎖状）
glc = 'C([C@H]([C@H]([C@@H]([C@H](C=O)O)O)O)O)O'
man = 'C([C@H]([C@H]([C@@H]([C@@H](C=O)O)O)O)O)O'
fru = 'C([C@H]([C@H]([C@@H](C(=O)CO)O)O)O)O'
red = lambda s: s.replace('C=O', 'CO')
glucitol, mannitol = can(red(glc)), can(red(man))
assert glucitol != mannitol
prods = set()
for chir in ('@', '@@'):
    s = fru.replace('C(=O)CO', f'[C{chir}H](O)CO')
    prods.add(can(s))
assert prods == {glucitol, mannitol}, prods
print('不斉炭素の数・還元生成物 OK', prods)

print('ALL OK')
