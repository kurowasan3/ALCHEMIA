# 検証スクリプト ORG-STR-023  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-STR-023/item.json'
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
    p = it['problem']
    n_items = len(re.findall(r'\\item\b', p.split('\\end{enumerate}')[0].split('\\begin{enumerate}')[-1])) if '\\begin{enumerate}' in p else 0
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert p.count('\\begin{enumerate}') <= 1
    assert n_items == len(nums) or (n_items == 0 and len(nums) == 0), (n_items, nums)
    assert nums == [str(i+1) for i in range(len(nums))], nums
    for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
        assert k in it, k
    print('小問数', n_items)
check_structure()
check_all_equations()
# ---- 個別の検証 ORG-STR-023：RDKit で立体異性体を列挙 ----
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
def info(smi):
    m = Chem.MolFromSmiles(smi)
    chir = Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=True)
    isos = {Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m, options=StereoEnumerationOptions(unique=True, onlyUnassigned=True))}
    def mirror(s):
        return Chem.CanonSmiles(s.replace('@@','!').replace('@','@@').replace('!','@'))
    meso = [s for s in isos if len(chir) and mirror(s) == Chem.CanonSmiles(s)]
    return len(chir), len(isos), meso
D = {'a':'ClC1CCCCC1Cl', 'b':'ClC1CCCC(Cl)C1', 'c':'ClC1CCC(Cl)CC1', 'd':'BrC1CCCCC1Cl'}
exp = {'a':(2,3,1), 'b':(2,3,1), 'c':(0,2,0), 'd':(2,4,0)}
for k, s in D.items():
    nc, ni, meso = info(s); print(k, s, '不斉炭素', nc, '立体異性体', ni, 'メソ', meso)
    assert (nc, ni, len(meso)) == exp[k], k
# (c) の 2 種がともに鏡像と同一（光学不活性）
m = Chem.MolFromSmiles(D['c'])
isos = [Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m)]
print('(c)', isos)
for s in isos:
    assert Chem.CanonSmiles(s.replace('@@','!').replace('@','@@').replace('!','@')) == Chem.CanonSmiles(s)
print('ALL OK')
