# 検証スクリプト ORG-HYD-025  python3 -I この.py
ITEM = '/home/claude/chem3/out/items/ORG-HYD-025/item.json'
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
# ---- 個別の検証 ORG-HYD-025 ----
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
def dou(c, h): return (2*c + 2 - h)//2
def rings(s): return Chem.MolFromSmiles(s).GetRingInfo().NumRings()
def pi(s):
    m = Chem.MolFromSmiles(s); n = 0
    for b in m.GetBonds():
        n += {Chem.BondType.DOUBLE:1, Chem.BondType.TRIPLE:2, Chem.BondType.AROMATIC:0.5}.get(b.GetBondType(), 0)
    return n
ex = {'A':('C1CCCCC1','C6H12',3), 'B':('C=CCCCC','C6H12',2), 'C':('C1=CCCCC1','C6H10',4),
      'D1':('C#CCCCC','C6H10',2), 'D2':('C=CCCC=C','C6H10',2), 'E':('Cc1ccccc1','C7H8',5), 'F':('C1CC2CC2C1','C6H10',3)}
for k,(s,f,cls) in ex.items():
    m = Chem.MolFromSmiles(s); assert CalcMolFormula(m) == f, (k, CalcMolFormula(m))
    c = int(f[1]); h = int(f.split('H')[1])
    assert rings(s) + pi(s) == dou(c,h), k
    # 水素付加量（芳香環は 3H2，C=C は 1，C≡C は 2）
    print(k, s, f, '不飽和度', dou(c,h), '環', rings(s), 'π', pi(s), '分類', cls)
# C: 付加 H2 1 mol → C6H12 環1（臭素水脱色しない）
assert dou(6,12) == 1 and rings('C1=CCCCC1') == 1
# D: 2 mol 付加で C6H14（不飽和度0）→ 環なし
assert dou(6,10) - 2 == dou(6,14) == 0
# E: 3 mol 付加で C7H14 不飽和度 1（環が残る）
assert dou(7,8) == 4 and dou(7,14) == 1
# F: 不飽和度 2 すべて環
assert rings('C1CC2CC2C1') == 2 and pi('C1CC2CC2C1') == 0
assert check_equation('C6H5CH3 + 3H2 -> C6H11CH3') is None
print('ALL OK')
