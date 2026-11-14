# 検証スクリプト ORG-ARO-044  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-ARO-044/item.json'
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
# ---- 個別の検証 ORG-ARO-044 ----
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import MolWt
W = {'C':12, 'H':1.0, 'O':16}
def M(f):
    import re
    return sum(W[e]*int(n or 1) for e, n in re.findall(r'([A-Z][a-z]?)(\d*)', f))
pa, ph, pp = 'O=C1OC(=O)c2ccccc21', 'Oc1ccccc1', 'O=C1OC(c2ccc(O)cc2)(c2ccc(O)cc2)c2ccccc21'
assert CalcMolFormula(Chem.MolFromSmiles(pa)) == 'C8H4O3' and CalcMolFormula(Chem.MolFromSmiles(pp)) == 'C20H14O4'
assert M('C8H4O3') == 148 and M('C6H6O') == 94 and M('C20H14O4') == 318
nA, nP = 0.296/148, 0.282/94
assert abs(nA-2.00e-3) < 1e-9 and abs(nP-3.00e-3) < 1e-9
prod = min(nA, nP/2); m = prod*318; print(prod, m); assert sig(m,3) == 0.477
left = nA - prod; assert sig(left,3) == 5.00e-4
assert sig(2*nA*94, 3) == 0.376
# (5) 反応位置：OH のパラ位に H があるか
for name, smi, exp in (('o-クレゾール', 'Cc1ccccc1O', True), ('p-クレゾール', 'Cc1ccc(O)cc1', False)):
    m_ = Chem.MolFromSmiles(smi)
    patt = Chem.MolFromSmarts('[OX2H]c1cc[cH]cc1')
    ok = m_.HasSubstructMatch(patt); print(name, 'パラ位に H', ok); assert ok == exp
print('ALL OK')
